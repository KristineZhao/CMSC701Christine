package main

import (
	"bytes"
	"compress/flate"
	"math"
	"sort"
)

// DeltaEncoder handles delta encoding between similar cells
type DeltaEncoder struct {
	threshold   float64
	quantLevels uint32
	lossy       bool
}

// NewDeltaEncoder creates a new delta encoder
func NewDeltaEncoder(lossy bool, threshold float64, quantLevels uint32) *DeltaEncoder {
	return &DeltaEncoder{
		threshold:   threshold,
		quantLevels: quantLevels,
		lossy:       lossy,
	}
}

// JaccardSimilarity calculates the Jaccard similarity between two gene sets
func JaccardSimilarity(genesA, genesB []uint32) float64 {
	if len(genesA) == 0 && len(genesB) == 0 {
		return 1.0
	}
	if len(genesA) == 0 || len(genesB) == 0 {
		return 0.0
	}

	// Create sets for faster lookup
	setA := make(map[uint32]bool)
	for _, gene := range genesA {
		setA[gene] = true
	}

	intersection := 0
	setB := make(map[uint32]bool)
	for _, gene := range genesB {
		if setA[gene] {
			intersection++
		}
		setB[gene] = true
	}

	union := len(setA) + len(setB) - intersection
	if union == 0 {
		return 1.0
	}

	return float64(intersection) / float64(union)
}

// FindBestReference finds the most similar cell to use as reference for delta encoding
func (de *DeltaEncoder) FindBestReference(targetCell SparseRow, candidates []SparseRow, candidateIndices []int) int {
	if len(candidates) == 0 {
		return -1
	}

	bestSimilarity := -1.0
	bestIndex := -1

	for i, candidate := range candidates {
		similarity := JaccardSimilarity(targetCell.Indices, candidate.Indices)
		if similarity > bestSimilarity && similarity > 0.1 { // Minimum similarity threshold
			bestSimilarity = similarity
			bestIndex = candidateIndices[i]
		}
	}

	return bestIndex
}

// ComputeDelta computes the delta between two sparse rows
func (de *DeltaEncoder) ComputeDelta(target, reference SparseRow) []int32 {
	// Create maps for faster lookup
	refMap := make(map[uint32]uint32)
	for i, gene := range reference.Indices {
		refMap[gene] = reference.Values[i]
	}

	targetMap := make(map[uint32]uint32)
	for i, gene := range target.Indices {
		targetMap[gene] = target.Values[i]
	}

	// Find all genes that appear in either target or reference
	allGenes := make(map[uint32]bool)
	for _, gene := range target.Indices {
		allGenes[gene] = true
	}
	for _, gene := range reference.Indices {
		allGenes[gene] = true
	}

	// Convert to sorted slice
	genes := make([]uint32, 0, len(allGenes))
	for gene := range allGenes {
		genes = append(genes, gene)
	}
	sort.Slice(genes, func(i, j int) bool { return genes[i] < genes[j] })

	// Compute deltas
	deltas := make([]int32, 0, len(genes))
	for _, gene := range genes {
		targetVal := targetMap[gene]
		refVal := refMap[gene]
		
		delta := int32(targetVal) - int32(refVal)
		
		// Apply lossy compression if enabled
		if de.lossy && math.Abs(float64(delta)) < de.threshold {
			delta = 0
		}
		
		deltas = append(deltas, delta)
	}

	return deltas
}

// QuantizeValue applies logarithmic quantization to a value
func (de *DeltaEncoder) QuantizeValue(value uint32) uint32 {
	if !de.lossy || value == 0 {
		return value
	}

	// Logarithmic quantization
	logVal := math.Log2(float64(value + 1))
	maxLog := math.Log2(float64(de.quantLevels))
	
	quantized := uint32(math.Round(logVal / maxLog * float64(de.quantLevels-1)))
	if quantized >= de.quantLevels {
		quantized = de.quantLevels - 1
	}

	// Convert back to approximate original scale
	return uint32(math.Pow(2, float64(quantized)*maxLog/float64(de.quantLevels-1))) - 1
}

// DequantizeValue reverses the quantization process
func (de *DeltaEncoder) DequantizeValue(quantized uint32) uint32 {
	if !de.lossy || quantized == 0 {
		return quantized
	}

	maxLog := math.Log2(float64(de.quantLevels))
	logVal := float64(quantized) * maxLog / float64(de.quantLevels-1)
	
	return uint32(math.Pow(2, logVal)) - 1
}

// CompressDeltas compresses a delta array using entropy coding
func (de *DeltaEncoder) CompressDeltas(deltas []int32) ([]byte, error) {
	if len(deltas) == 0 {
		return []byte{}, nil
	}

	var buf bytes.Buffer
	
	// Use flate compression (DEFLATE algorithm)
	writer, err := flate.NewWriter(&buf, flate.BestCompression)
	if err != nil {
		return nil, err
	}

	// Convert int32 deltas to bytes using variable-length encoding
	for _, delta := range deltas {
		if err := writeVarint(&buf, delta); err != nil {
			writer.Close()
			return nil, err
		}
	}

	writer.Close()
	return buf.Bytes(), nil
}

// DecompressDeltas decompresses a delta array
func (de *DeltaEncoder) DecompressDeltas(compressed []byte) ([]int32, error) {
	if len(compressed) == 0 {
		return []int32{}, nil
	}

	buf := bytes.NewReader(compressed)
	reader := flate.NewReader(buf)
	defer reader.Close()

	var deltas []int32
	decompressedBuf := bytes.NewBuffer(nil)
	if _, err := decompressedBuf.ReadFrom(reader); err != nil {
		return nil, err
	}

	decompressedReader := bytes.NewReader(decompressedBuf.Bytes())
	for decompressedReader.Len() > 0 {
		delta, err := readVarint(decompressedReader)
		if err != nil {
			break // End of data
		}
		deltas = append(deltas, delta)
	}

	return deltas, nil
}

// ReconstructFromDelta reconstructs the target cell from reference and delta
func (de *DeltaEncoder) ReconstructFromDelta(reference SparseRow, deltas []int32, geneIndices []uint32) SparseRow {
	// Create reference map
	refMap := make(map[uint32]uint32)
	for i, gene := range reference.Indices {
		refMap[gene] = reference.Values[i]
	}

	// Apply deltas
	var resultIndices []uint32
	var resultValues []uint32

	for i, gene := range geneIndices {
		if i >= len(deltas) {
			break
		}
		
		refVal := refMap[gene]
		newVal := int32(refVal) + deltas[i]
		
		if newVal > 0 {
			resultIndices = append(resultIndices, gene)
			resultValues = append(resultValues, uint32(newVal))
		}
	}

	return SparseRow{
		Indices: resultIndices,
		Values:  resultValues,
	}
}

// writeVarint writes a signed integer using variable-length encoding
func writeVarint(buf *bytes.Buffer, value int32) error {
	// Zigzag encoding to handle signed integers
	uvalue := uint32((value << 1) ^ (value >> 31))
	
	for uvalue >= 0x80 {
		buf.WriteByte(byte(uvalue | 0x80))
		uvalue >>= 7
	}
	buf.WriteByte(byte(uvalue))
	return nil
}

// readVarint reads a variable-length encoded signed integer
func readVarint(buf *bytes.Reader) (int32, error) {
	var uvalue uint32
	var shift uint
	
	for {
		b, err := buf.ReadByte()
		if err != nil {
			return 0, err
		}
		
		uvalue |= uint32(b&0x7F) << shift
		if b < 0x80 {
			break
		}
		shift += 7
		if shift >= 32 {
			return 0, bytes.ErrTooLarge
		}
	}
	
	// Zigzag decoding
	value := int32((uvalue >> 1) ^ (-(uvalue & 1)))
	return value, nil
}
