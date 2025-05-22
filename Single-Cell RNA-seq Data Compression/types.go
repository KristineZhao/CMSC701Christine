package main

import (
	"encoding/binary"
	"io"
)

// SparseRow represents a single cell's expression profile
type SparseRow struct {
	Indices []uint32 // Gene indices (sorted)
	Values  []uint32 // Expression counts
}

// CompressedData represents the complete compressed dataset
type CompressedData struct {
	Header       Header
	GeneNames    []string
	CellNames    []string
	CompressedRows []CompressedRow
}

// Header contains metadata about the compressed data
type Header struct {
	Version      uint32
	NumCells     uint32
	NumGenes     uint32
	IsLossy      bool
	Threshold    float64
	QuantLevels  uint32
	Timestamp    int64
}

// CompressedRow represents a compressed cell's expression profile
type CompressedRow struct {
	EliasGenes   []byte  // Elias-Fano encoded gene indices
	DeltaValues  []byte  // Delta-encoded and compressed expression values
	RefCell      int32   // Reference cell index for delta encoding (-1 if none)
	NumGenes     uint32  // Number of expressed genes
	MaxGeneIndex uint32  // Maximum gene index for Elias-Fano
}

// EliasRange represents the range information for Elias-Fano encoding
type EliasRange struct {
	Universe uint32 // Universe size (u)
	Count    uint32 // Number of elements (k)
	LowBits  uint32 // Number of low bits (l)
}

// BitArray provides efficient bit-level operations
type BitArray struct {
	Data []uint64
	Size uint32
}

// NewBitArray creates a new bit array with the specified size
func NewBitArray(size uint32) *BitArray {
	numWords := (size + 63) / 64
	return &BitArray{
		Data: make([]uint64, numWords),
		Size: size,
	}
}

// SetBit sets the bit at the given position
func (ba *BitArray) SetBit(pos uint32) {
	if pos >= ba.Size {
		return
	}
	wordIndex := pos / 64
	bitIndex := pos % 64
	ba.Data[wordIndex] |= 1 << bitIndex
}

// GetBit returns the bit at the given position
func (ba *BitArray) GetBit(pos uint32) bool {
	if pos >= ba.Size {
		return false
	}
	wordIndex := pos / 64
	bitIndex := pos % 64
	return (ba.Data[wordIndex] & (1 << bitIndex)) != 0
}

// WriteBits writes the specified number of bits starting from the given position
func (ba *BitArray) WriteBits(pos uint32, value uint64, numBits uint32) {
	for i := uint32(0); i < numBits; i++ {
		if (value & (1 << i)) != 0 {
			ba.SetBit(pos + i)
		}
	}
}

// ReadBits reads the specified number of bits starting from the given position
func (ba *BitArray) ReadBits(pos uint32, numBits uint32) uint64 {
	var result uint64
	for i := uint32(0); i < numBits; i++ {
		if ba.GetBit(pos + i) {
			result |= 1 << i
		}
	}
	return result
}

// WriteTo writes the bit array to an io.Writer
func (ba *BitArray) WriteTo(w io.Writer) error {
	// Write size first
	if err := binary.Write(w, binary.LittleEndian, ba.Size); err != nil {
		return err
	}
	// Write data
	return binary.Write(w, binary.LittleEndian, ba.Data)
}

// ReadFrom reads the bit array from an io.Reader
func (ba *BitArray) ReadFrom(r io.Reader) error {
	// Read size
	if err := binary.Read(r, binary.LittleEndian, &ba.Size); err != nil {
		return err
	}
	// Calculate number of words needed
	numWords := (ba.Size + 63) / 64
	ba.Data = make([]uint64, numWords)
	// Read data
	return binary.Read(r, binary.LittleEndian, ba.Data)
}

// CellSimilarity represents similarity between two cells
type CellSimilarity struct {
	CellA      uint32
	CellB      uint32
	Similarity float64
}

// CompressionStats holds statistics about compression performance
type CompressionStats struct {
	OriginalSize    int64
	CompressedSize  int64
	CompressionTime int64 // nanoseconds
	CompressionRatio float64
	NumCells        uint32
	NumGenes        uint32
	AvgGenesPerCell float64
	Sparsity        float64
}

// DecompressionStats holds statistics about decompression performance
type DecompressionStats struct {
	DecompressedSize   int64
	DecompressionTime  int64 // nanoseconds
	NumCells           uint32
	NumGenes           uint32
	AvgGenesPerCell    float64
}
