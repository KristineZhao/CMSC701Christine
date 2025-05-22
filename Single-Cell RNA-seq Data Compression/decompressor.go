package main

import (
	"fmt"
	"runtime"
	"sync"
	"time"
)

// Decompressor handles the decompression of single-cell RNA-seq data
type Decompressor struct{}

// NewDecompressor creates a new decompressor
func NewDecompressor() *Decompressor {
	return &Decompressor{}
}

// Decompress decompresses the compressed data back to sparse matrix format
func (d *Decompressor) Decompress(compressed *CompressedData) ([]SparseRow, []string, []string, error) {
	startTime := time.Now()

	matrix := make([]SparseRow, compressed.Header.NumCells)
	numWorkers := runtime.NumCPU()
	jobs := make(chan int, compressed.Header.NumCells)
	var wg sync.WaitGroup
	var mu sync.Mutex
	var decompressErr error

	// Create delta encoder for decompression
	deltaEncoder := NewDeltaEncoder(
		compressed.Header.IsLossy,
		compressed.Header.Threshold,
		compressed.Header.QuantLevels,
	)

	// Start workers
	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for cellIdx := range jobs {
				row, err := d.decompressCell(
					compressed.CompressedRows[cellIdx],
					matrix,
					deltaEncoder,
					compressed.Header.NumGenes,
				)
				if err != nil {
					mu.Lock()
					if decompressErr == nil {
						decompressErr = fmt.Errorf("error decompressing cell %d: %w", cellIdx, err)
					}
					mu.Unlock()
					continue
				}
				
				mu.Lock()
				matrix[cellIdx] = row
				mu.Unlock()
			}
		}()
	}

	// Send jobs
	for i := 0; i < int(compressed.Header.NumCells); i++ {
		jobs <- i
	}
	close(jobs)
	wg.Wait()

	if decompressErr != nil {
		return nil, nil, nil, decompressErr
	}

	// Apply dequantization if lossy compression was used
	if compressed.Header.IsLossy {
		matrix = d.applyDequantization(matrix, deltaEncoder)
	}

	fmt.Printf("Decompression completed in %v\n", time.Since(startTime))
	return matrix, compressed.GeneNames, compressed.CellNames, nil
}

// decompressCell decompresses a single cell's expression profile
func (d *Decompressor) decompressCell(
	compressedRow CompressedRow,
	matrix []SparseRow,
	deltaEncoder *DeltaEncoder,
	numGenes uint32,
) (SparseRow, error) {
	var result SparseRow

	// Decompress gene indices using Elias-Fano decoding
	if len(compressedRow.EliasGenes) > 0 {
		decoder, err := NewEliasDecoder(compressedRow.EliasGenes)
		if err != nil {
			return result, fmt.Errorf("failed to create Elias-Fano decoder: %w", err)
		}

		geneIndices, err := decoder.Decode()
		if err != nil {
			return result, fmt.Errorf("failed to decode gene indices: %w", err)
		}

		result.Indices = geneIndices
	}

	// Decompress expression values
	if len(compressedRow.DeltaValues) > 0 {
		deltas, err := deltaEncoder.DecompressDeltas(compressedRow.DeltaValues)
		if err != nil {
			return result, fmt.Errorf("failed to decompress deltas: %w", err)
		}

		if compressedRow.RefCell >= 0 && int(compressedRow.RefCell) < len(matrix) {
			// Reconstruct using reference cell and deltas
			refCell := matrix[compressedRow.RefCell]
			result = deltaEncoder.ReconstructFromDelta(refCell, deltas, result.Indices)
		} else {
			// No reference cell, deltas are the actual values
			if len(deltas) != len(result.Indices) {
				return result, fmt.Errorf("mismatch between number of genes (%d) and values (%d)", 
					len(result.Indices), len(deltas))
			}

			result.Values = make([]uint32, len(deltas))
			for i, delta := range deltas {
				if delta < 0 {
					result.Values[i] = 0
				} else {
					result.Values[i] = uint32(delta)
				}
			}
		}
	}

	// Ensure indices and values have the same length
	if len(result.Indices) != len(result.Values) {
		minLen := len(result.Indices)
		if len(result.Values) < minLen {
			minLen = len(result.Values)
		}
		result.Indices = result.Indices[:minLen]
		result.Values = result.Values[:minLen]
	}

	return result, nil
}

// applyDequantization applies dequantization to restore approximate original values
func (d *Decompressor) applyDequantization(matrix []SparseRow, deltaEncoder *DeltaEncoder) []SparseRow {
	if !deltaEncoder.lossy {
		return matrix
	}

	dequantized := make([]SparseRow, len(matrix))
	for i, row := range matrix {
		dequantizedValues := make([]uint32, len(row.Values))
		for j, value := range row.Values {
			dequantizedValues[j] = deltaEncoder.DequantizeValue(value)
		}
		dequantized[i] = SparseRow{
			Indices: append([]uint32(nil), row.Indices...),
			Values:  dequantizedValues,
		}
	}
	return dequantized
}
