package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"compress/zlib"
	"encoding/binary"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// LoadSparseMatrix loads a sparse matrix from various file formats
func LoadSparseMatrix(filename string) ([]SparseRow, []string, []string, error) {
	ext := strings.ToLower(filename[strings.LastIndex(filename, "."):])
	
	switch ext {
	case ".csv", ".tsv":
		return loadFromCSV(filename, ext == ".tsv")
	case ".gz":
		// Handle compressed files
		if strings.HasSuffix(strings.ToLower(filename), ".csv.gz") {
			return loadFromCompressedCSV(filename, false)
		} else if strings.HasSuffix(strings.ToLower(filename), ".tsv.gz") {
			return loadFromCompressedCSV(filename, true)
		}
		return nil, nil, nil, fmt.Errorf("unsupported compressed file format: %s", filename)
	case ".rds":
		return loadFromRDS(filename)
	default:
		return nil, nil, nil, fmt.Errorf("unsupported file format: %s", ext)
	}
}

// loadFromCSV loads matrix data from CSV/TSV files
func loadFromCSV(filename string, isTab bool) ([]SparseRow, []string, []string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, nil, err
	}
	defer file.Close()

	return parseCSVReader(file, isTab)
}

// loadFromCompressedCSV loads matrix data from compressed CSV/TSV files
func loadFromCompressedCSV(filename string, isTab bool) ([]SparseRow, []string, []string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, nil, err
	}
	defer file.Close()

	gzReader, err := gzip.NewReader(file)
	if err != nil {
		return nil, nil, nil, err
	}
	defer gzReader.Close()

	return parseCSVReader(gzReader, isTab)
}

// parseCSVReader parses CSV data from an io.Reader
func parseCSVReader(reader io.Reader, isTab bool) ([]SparseRow, []string, []string, error) {
	csvReader := csv.NewReader(reader)
	if isTab {
		csvReader.Comma = '\t'
	}

	// Read header (gene names)
	header, err := csvReader.Read()
	if err != nil {
		return nil, nil, nil, fmt.Errorf("failed to read header: %w", err)
	}

	// First column is usually cell names, rest are gene names
	geneNames := header[1:]
	var cellNames []string
	var matrix []SparseRow

	// Read data rows
	for {
		record, err := csvReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, nil, fmt.Errorf("failed to read CSV record: %w", err)
		}

		if len(record) < 2 {
			continue // Skip invalid rows
		}

		cellName := record[0]
		cellNames = append(cellNames, cellName)

		// Parse expression values
		var indices []uint32
		var values []uint32

		for i, valueStr := range record[1:] {
			if valueStr == "" || valueStr == "0" {
				continue // Skip zero values
			}

			value, err := strconv.ParseFloat(valueStr, 64)
			if err != nil {
				continue // Skip invalid values
			}

			if value > 0 {
				indices = append(indices, uint32(i))
				values = append(values, uint32(value))
			}
		}

		matrix = append(matrix, SparseRow{
			Indices: indices,
			Values:  values,
		})
	}

	return matrix, geneNames, cellNames, nil
}

// loadFromRDS loads matrix data from RDS files (simplified implementation)
// Note: This is a basic implementation and may not handle all RDS formats
func loadFromRDS(filename string) ([]SparseRow, []string, []string, error) {
	// For now, return an error suggesting conversion to CSV
	return nil, nil, nil, fmt.Errorf("RDS format not fully supported yet. Please convert to CSV/TSV format using R:\n" +
		"library(Matrix)\n" +
		"data <- readRDS('%s')\n" +
		"write.csv(as.matrix(data$all_data[[1]]$hg19$mat), 'output.csv')", filename)
}

// SaveSparseMatrix saves a sparse matrix to a CSV file
func SaveSparseMatrix(matrix []SparseRow, geneNames, cellNames []string, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write header
	header := append([]string{"Cell"}, geneNames...)
	if err := writer.Write(header); err != nil {
		return err
	}

	// Write data rows
	for i, row := range matrix {
		cellName := ""
		if i < len(cellNames) {
			cellName = cellNames[i]
		} else {
			cellName = fmt.Sprintf("Cell_%d", i+1)
		}

		// Create dense row with zeros
		denseRow := make([]string, len(geneNames)+1)
		denseRow[0] = cellName
		for j := 1; j < len(denseRow); j++ {
			denseRow[j] = "0"
		}

		// Fill in non-zero values
		for j, geneIdx := range row.Indices {
			if int(geneIdx) < len(geneNames) {
				denseRow[geneIdx+1] = strconv.Itoa(int(row.Values[j]))
			}
		}

		if err := writer.Write(denseRow); err != nil {
			return err
		}
	}

	return nil
}

// SaveToFile saves compressed data to a binary file
func (cd *CompressedData) SaveToFile(filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	// Use zlib compression for the entire file
	zlibWriter := zlib.NewWriter(file)
	defer zlibWriter.Close()

	var buf bytes.Buffer

	// Write header
	if err := binary.Write(&buf, binary.LittleEndian, cd.Header); err != nil {
		return err
	}

	// Write gene names
	if err := writeStringSlice(&buf, cd.GeneNames); err != nil {
		return err
	}

	// Write cell names
	if err := writeStringSlice(&buf, cd.CellNames); err != nil {
		return err
	}

	// Write number of compressed rows
	if err := binary.Write(&buf, binary.LittleEndian, uint32(len(cd.CompressedRows))); err != nil {
		return err
	}

	// Write compressed rows
	for _, row := range cd.CompressedRows {
		if err := writeCompressedRow(&buf, row); err != nil {
			return err
		}
	}

	// Write buffer to zlib writer
	_, err = zlibWriter.Write(buf.Bytes())
	return err
}

// LoadCompressedData loads compressed data from a binary file
func LoadCompressedData(filename string) (*CompressedData, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	// Use zlib decompression
	zlibReader, err := zlib.NewReader(file)
	if err != nil {
		return nil, err
	}
	defer zlibReader.Close()

	var buf bytes.Buffer
	if _, err := buf.ReadFrom(zlibReader); err != nil {
		return nil, err
	}

	reader := bytes.NewReader(buf.Bytes())
	cd := &CompressedData{}

	// Read header
	if err := binary.Read(reader, binary.LittleEndian, &cd.Header); err != nil {
		return nil, err
	}

	// Read gene names
	cd.GeneNames, err = readStringSlice(reader)
	if err != nil {
		return nil, err
	}

	// Read cell names
	cd.CellNames, err = readStringSlice(reader)
	if err != nil {
		return nil, err
	}

	// Read number of compressed rows
	var numRows uint32
	if err := binary.Read(reader, binary.LittleEndian, &numRows); err != nil {
		return nil, err
	}

	// Read compressed rows
	cd.CompressedRows = make([]CompressedRow, numRows)
	for i := uint32(0); i < numRows; i++ {
		row, err := readCompressedRow(reader)
		if err != nil {
			return nil, err
		}
		cd.CompressedRows[i] = row
	}

	return cd, nil
}

// EstimateSize estimates the size of compressed data in bytes
func (cd *CompressedData) EstimateSize() int {
	size := 0
	
	// Header size
	size += 32 // approximate header size
	
	// Gene names
	for _, name := range cd.GeneNames {
		size += len(name) + 4 // string length + length prefix
	}
	
	// Cell names
	for _, name := range cd.CellNames {
		size += len(name) + 4 // string length + length prefix
	}
	
	// Compressed rows
	for _, row := range cd.CompressedRows {
		size += len(row.EliasGenes) + len(row.DeltaValues) + 16 // data + metadata
	}
	
	return size
}

// Helper functions for reading/writing binary data

func writeStringSlice(buf *bytes.Buffer, strings []string) error {
	// Write number of strings
	if err := binary.Write(buf, binary.LittleEndian, uint32(len(strings))); err != nil {
		return err
	}
	
	// Write each string
	for _, s := range strings {
		if err := writeString(buf, s); err != nil {
			return err
		}
	}
	return nil
}

func readStringSlice(reader *bytes.Reader) ([]string, error) {
	var count uint32
	if err := binary.Read(reader, binary.LittleEndian, &count); err != nil {
		return nil, err
	}
	
	strings := make([]string, count)
	for i := uint32(0); i < count; i++ {
		s, err := readString(reader)
		if err != nil {
			return nil, err
		}
		strings[i] = s
	}
	return strings, nil
}

func writeString(buf *bytes.Buffer, s string) error {
	// Write string length
	if err := binary.Write(buf, binary.LittleEndian, uint32(len(s))); err != nil {
		return err
	}
	// Write string data
	_, err := buf.WriteString(s)
	return err
}

func readString(reader *bytes.Reader) (string, error) {
	var length uint32
	if err := binary.Read(reader, binary.LittleEndian, &length); err != nil {
		return "", err
	}
	
	data := make([]byte, length)
	if _, err := reader.Read(data); err != nil {
		return "", err
	}
	return string(data), nil
}

func writeCompressedRow(buf *bytes.Buffer, row CompressedRow) error {
	// Write metadata
	if err := binary.Write(buf, binary.LittleEndian, row.RefCell); err != nil {
		return err
	}
	if err := binary.Write(buf, binary.LittleEndian, row.NumGenes); err != nil {
		return err
	}
	if err := binary.Write(buf, binary.LittleEndian, row.MaxGeneIndex); err != nil {
		return err
	}
	
	// Write Elias-Fano data
	if err := binary.Write(buf, binary.LittleEndian, uint32(len(row.EliasGenes))); err != nil {
		return err
	}
	if _, err := buf.Write(row.EliasGenes); err != nil {
		return err
	}
	
	// Write delta values
	if err := binary.Write(buf, binary.LittleEndian, uint32(len(row.DeltaValues))); err != nil {
		return err
	}
	_, err := buf.Write(row.DeltaValues)
	return err
}

func readCompressedRow(reader *bytes.Reader) (CompressedRow, error) {
	var row CompressedRow
	
	// Read metadata
	if err := binary.Read(reader, binary.LittleEndian, &row.RefCell); err != nil {
		return row, err
	}
	if err := binary.Read(reader, binary.LittleEndian, &row.NumGenes); err != nil {
		return row, err
	}
	if err := binary.Read(reader, binary.LittleEndian, &row.MaxGeneIndex); err != nil {
		return row, err
	}
	
	// Read Elias-Fano data
	var eliasLen uint32
	if err := binary.Read(reader, binary.LittleEndian, &eliasLen); err != nil {
		return row, err
	}
	row.EliasGenes = make([]byte, eliasLen)
	if _, err := reader.Read(row.EliasGenes); err != nil {
		return row, err
	}
	
	// Read delta values
	var deltaLen uint32
	if err := binary.Read(reader, binary.LittleEndian, &deltaLen); err != nil {
		return row, err
	}
	row.DeltaValues = make([]byte, deltaLen)
	if _, err := reader.Read(row.DeltaValues); err != nil {
		return row, err
	}
	
	return row, nil
}
