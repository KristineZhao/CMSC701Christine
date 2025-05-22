package main

import (
	"fmt"
	"runtime"
	"sort"
	"sync"
	"time"
)

// Compressor handles the compression of single-cell RNA-seq data
type Compressor struct {
	lossy       bool
	threshold   float64
	quantLevels uint32
	deltaEncoder *DeltaEncoder
}

// NewCompressor creates a new compressor with the specified parameters
func NewCompressor(lossy bool, threshold float64, quantLevels uint32) *Compressor {
	return &Compressor{
		lossy:       lossy,
		threshold:   threshold,
		quantLevels: quantLevels,
		deltaEncoder: NewDeltaEncoder(lossy, threshold, quantLevels),
	}
