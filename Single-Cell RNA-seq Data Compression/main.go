package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"
)

func main() {
	var (
		inputFile    = flag.String("input", "", "Input file path (CSV, TSV, or RDS)")
		outputFile   = flag.String("output", "", "Output compressed file path")
		mode         = flag.String("mode", "compress", "Mode: compress or decompress")
		lossy        = flag.Bool("lossy", false, "Enable lossy compression")
		threshold    = flag.Float64("threshold", 0.1, "Delta threshold for lossy compression")
		quantLevels  = flag.Int("quant", 256, "Quantization levels for lossy compression")
		verbose      = flag.Bool("verbose", false, "Verbose output")
	)
	flag.Parse()

	if *inputFile == "" {
		fmt.Println("Usage:")
		fmt.Println("  Compress: go run . -input data.csv -output compressed.scz -mode compress")
		fmt.Println("  Decompress: go run . -input compressed.scz -output decompressed.csv -mode decompress")
		fmt.Println("  Lossy: go run . -input data.csv -output compressed.scz -lossy -threshold 0.1")
		os.Exit(1)
	}

	switch *mode {
	case "compress":
		if *outputFile == "" {
			*outputFile = strings.TrimSuffix(*inputFile, filepath.Ext(*inputFile)) + ".scz"
		}
		err := compressFile(*inputFile, *outputFile, *lossy, *threshold, *quantLevels, *verbose)
		if err != nil {
			log.Fatalf("Compression failed: %v", err)
		}
		fmt.Printf("Successfully compressed %s to %s\n", *inputFile, *outputFile)

	case "decompress":
		if *outputFile == "" {
			*outputFile = strings.TrimSuffix(*inputFile, filepath.Ext(*inputFile)) + "_decompressed.csv"
		}
		err := decompressFile(*inputFile, *outputFile, *verbose)
		if err != nil {
			log.Fatalf("Decompression failed: %v", err)
		}
		fmt.Printf("Successfully decompressed %s to %s\n", *inputFile, *outputFile)

	default:
		log.Fatalf("Unknown mode: %s. Use 'compress' or 'decompress'", *mode)
	}
}

func compressFile(inputFile, outputFile string, lossy bool, threshold float64, quantLevels int, verbose bool) error {
	// Load the sparse matrix
	matrix, geneNames, cellNames, err := LoadSparseMatrix(inputFile)
	if err != nil {
		return fmt.Errorf("failed to load input file: %w", err)
	}

	if verbose {
		fmt.Printf("Loaded matrix: %d cells x %d genes\n", len(matrix), len(geneNames))
		fmt.Printf("Total non-zero entries: %d\n", countNonZeros(matrix))
	}

	// Create compressor
	compressor := NewCompressor(lossy, threshold, quantLevels)

	// Compress the matrix
	compressed, err := compressor.Compress(matrix, geneNames, cellNames)
	if err != nil {
		return fmt.Errorf("compression failed: %w", err)
	}

	// Save compressed data
	err = compressed.SaveToFile(outputFile)
	if err != nil {
		return fmt.Errorf("failed to save compressed file: %w", err)
	}

	if verbose {
		originalSize := estimateOriginalSize(matrix, geneNames, cellNames)
		compressedSize := compressed.EstimateSize()
		ratio := float64(originalSize) / float64(compressedSize)
		fmt.Printf("Original size: %d bytes\n", originalSize)
		fmt.Printf("Compressed size: %d bytes\n", compressedSize)
		fmt.Printf("Compression ratio: %.2fx\n", ratio)
	}

	return nil
}

func decompressFile(inputFile, outputFile string, verbose bool) error {
	// Load compressed data
	compressed, err := LoadCompressedData(inputFile)
	if err != nil {
		return fmt.Errorf("failed to load compressed file: %w", err)
	}

	// Create decompressor
	decompressor := NewDecompressor()

	// Decompress the data
	matrix, geneNames, cellNames, err := decompressor.Decompress(compressed)
	if err != nil {
		return fmt.Errorf("decompression failed: %w", err)
	}

	if verbose {
		fmt.Printf("Decompressed matrix: %d cells x %d genes\n", len(matrix), len(geneNames))
		fmt.Printf("Total non-zero entries: %d\n", countNonZeros(matrix))
	}

	// Save decompressed matrix
	err = SaveSparseMatrix(matrix, geneNames, cellNames, outputFile)
	if err != nil {
		return fmt.Errorf("failed to save decompressed file: %w", err)
	}

	return nil
}

func countNonZeros(matrix []SparseRow) int {
	count := 0
	for _, row := range matrix {
		count += len(row.Values)
	}
	return count
}

func estimateOriginalSize(matrix []SparseRow, geneNames, cellNames []string) int {
	size := 0
	// Gene names
	for _, name := range geneNames {
		size += len(name)
	}
	// Cell names
	for _, name := range cellNames {
		size += len(name)
	}
	// Matrix data (assuming 4 bytes per int32)
	for _, row := range matrix {
		size += len(row.Indices) * 4  // gene indices
		size += len(row.Values) * 4   // expression values
	}
	return size
}
