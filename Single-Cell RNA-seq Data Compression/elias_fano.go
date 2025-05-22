package main

import (
	"bytes"
	"encoding/binary"
	"fmt"
)

// EliasEncoder handles Elias-Fano encoding of sorted integer sequences
type EliasEncoder struct {
	universe uint32
	count    uint32
	lowBits  uint32
}

// NewEliasEncoder creates a new Elias-Fano encoder
func NewEliasEncoder(universe, count uint32) *EliasEncoder {
	lowBits := uint32(0)
	if count > 0 && universe > count {
		// Calculate l = floor(log2(u/k))
		ratio := universe / count
		for (1 << lowBits) < ratio {
			lowBits++
		}
		if lowBits > 0 {
			lowBits--
		}
	}

	return &EliasEncoder{
		universe: universe,
		count:    count,
		lowBits:  lowBits,
	}
}

// Encode compresses a sorted sequence of integers using Elias-Fano encoding
func (e *EliasEncoder) Encode(sequence []uint32) ([]byte, error) {
	if len(sequence) == 0 {
		return []byte{}, nil
	}

	if uint32(len(sequence)) != e.count {
		return nil, fmt.Errorf("sequence length %d doesn't match expected count %d", len(sequence), e.count)
	}

	// Validate that sequence is sorted and within universe
	for i, val := range sequence {
		if val >= e.universe {
			return nil, fmt.Errorf("value %d at index %d exceeds universe %d", val, i, e.universe)
		}
		if i > 0 && val <= sequence[i-1] {
			return nil, fmt.Errorf("sequence not sorted at index %d: %d <= %d", i, val, sequence[i-1])
		}
	}

	var buf bytes.Buffer

	// Write header information
	binary.Write(&buf, binary.LittleEndian, e.universe)
	binary.Write(&buf, binary.LittleEndian, e.count)
	binary.Write(&buf, binary.LittleEndian, e.lowBits)

	if e.count == 0 {
		return buf.Bytes(), nil
	}

	// Encode low bits
	lowBitsSize := e.count * e.lowBits
	lowArray := NewBitArray(lowBitsSize)
	
	for i, val := range sequence {
		lowValue := val & ((1 << e.lowBits) - 1)
		lowArray.WriteBits(uint32(i)*e.lowBits, uint64(lowValue), e.lowBits)
	}

	// Encode high bits using unary encoding
	highBitsSize := e.count + (e.universe >> e.lowBits) + 1
	highArray := NewBitArray(highBitsSize)
	
	pos := uint32(0)
	for _, val := range sequence {
		highValue := val >> e.lowBits
		pos += highValue
		if pos < highBitsSize {
			highArray.SetBit(pos)
		}
		pos++
	}

	// Write low bits array
	if err := lowArray.WriteTo(&buf); err != nil {
		return nil, err
	}

	// Write high bits array
	if err := highArray.WriteTo(&buf); err != nil {
		return nil, err
	}

	return buf.Bytes(), nil
}

// EliasDecoder handles Elias-Fano decoding
type EliasDecoder struct {
	universe uint32
	count    uint32
	lowBits  uint32
	lowArray *BitArray
	highArray *BitArray
}

// NewEliasDecoder creates a new Elias-Fano decoder from encoded data
func NewEliasDecoder(data []byte) (*EliasDecoder, error) {
	if len(data) < 12 { // 3 uint32s
		return nil, fmt.Errorf("encoded data too short")
	}

	buf := bytes.NewReader(data)
	decoder := &EliasDecoder{}

	// Read header
	if err := binary.Read(buf, binary.LittleEndian, &decoder.universe); err != nil {
		return nil, err
	}
	if err := binary.Read(buf, binary.LittleEndian, &decoder.count); err != nil {
		return nil, err
	}
	if err := binary.Read(buf, binary.LittleEndian, &decoder.lowBits); err != nil {
		return nil, err
	}

	if decoder.count == 0 {
		return decoder, nil
	}

	// Read low bits array
	decoder.lowArray = &BitArray{}
	if err := decoder.lowArray.ReadFrom(buf); err != nil {
		return nil, err
	}

	// Read high bits array
	decoder.highArray = &BitArray{}
	if err := decoder.highArray.ReadFrom(buf); err != nil {
		return nil, err
	}

	return decoder, nil
}

// Decode decompresses the Elias-Fano encoded sequence
func (d *EliasDecoder) Decode() ([]uint32, error) {
	if d.count == 0 {
		return []uint32{}, nil
	}

	result := make([]uint32, d.count)
	
	// Decode high bits to find positions
	highPos := uint32(0)
	currentHigh := uint32(0)
	
	for i := uint32(0); i < d.count; i++ {
		// Find next set bit in high array
		for highPos < d.highArray.Size && !d.highArray.GetBit(highPos) {
			highPos++
			currentHigh++
		}
		
		if highPos >= d.highArray.Size {
			return nil, fmt.Errorf("unexpected end of high bits array")
		}
		
		// Get low bits for this element
		lowValue := d.lowArray.ReadBits(i*d.lowBits, d.lowBits)
		
		// Combine high and low parts
		result[i] = (currentHigh << d.lowBits) | uint32(lowValue)
		
		highPos++
	}

	return result, nil
}

// Access provides random access to the i-th element without full decoding
func (d *EliasDecoder) Access(index uint32) (uint32, error) {
	if index >= d.count {
		return 0, fmt.Errorf("index %d out of range [0, %d)", index, d.count)
	}

	if d.count == 0 {
		return 0, fmt.Errorf("empty sequence")
	}

	// Get low bits for this element
	lowValue := d.lowArray.ReadBits(index*d.lowBits, d.lowBits)

	// Find the corresponding high value by counting set bits
	highValue := d.findHighValue(index)

	return (highValue << d.lowBits) | uint32(lowValue), nil
}

// findHighValue finds the high value for the given index by scanning the high bits array
func (d *EliasDecoder) findHighValue(index uint32) uint32 {
	setBitsCount := uint32(0)
	pos := uint32(0)
	
	for pos < d.highArray.Size && setBitsCount <= index {
		if d.highArray.GetBit(pos) {
			if setBitsCount == index {
				// Count the number of 0s before this position
				zeros := pos - setBitsCount
				return zeros
			}
			setBitsCount++
		}
		pos++
	}
	
	return 0 // Should not reach here in valid data
}

// Size returns the number of elements in the encoded sequence
func (d *EliasDecoder) Size() uint32 {
	return d.count
}

// Universe returns the universe size of the encoded sequence
func (d *EliasDecoder) Universe() uint32 {
	return d.universe
}
