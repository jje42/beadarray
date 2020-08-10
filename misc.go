package beadarray

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
)

func readNextBytes(file io.Reader, number int) ([]byte, error) {
	buf := make([]byte, number)
	n, err := io.ReadAtLeast(file, buf, number)
	if err != nil {
		return buf, err
	}
	if n != number {
		return buf, fmt.Errorf("readbytes expected %v got %v", number, n)
	}
	return buf, nil
}

func readInt(file io.Reader) (int, error) {
	b, err := readNextBytes(file, 4)
	if err != nil {
		return 0, fmt.Errorf("readInt failed reading bytes: %w", err)
	}
	var r int32
	buf := bytes.NewReader(b)
	err = binary.Read(buf, binary.LittleEndian, &r)
	if err != nil {
		return 0, fmt.Errorf("readInt failed: %w", err)
	}
	return int(r), nil
}

func mustReadInt(r io.Reader) int {
	n, err := readInt(r)
	if err != nil {
		panic(err)
	}
	return n
}

func readUint8(file io.Reader) (uint8, error) {
	var x uint8
	b, err := readNextBytes(file, 1)
	if err != nil {
		return 0, err
	}
	buf := bytes.NewReader(b)
	err = binary.Read(buf, binary.LittleEndian, &x)
	if err != nil {
		return 0, err
	}
	return x, nil
}

func readInt16(file io.Reader) (int16, error) {
	var r int16
	// b, err := readNextBytes(file, 2)
	// if err != nil {
	// 	return r, err
	// }
	// buf := bytes.NewReader(b)
	// err = binary.Read(buf, binary.LittleEndian, &r)
	err := binary.Read(file, binary.LittleEndian, &r)
	return r, err
}

func readUint16(r io.Reader) (uint16, error) {
	var ret uint16
	b, err := readNextBytes(r, 2)
	if err != nil {
		return ret, err
	}
	buf := bytes.NewReader(b)
	err = binary.Read(buf, binary.LittleEndian, &r)
	return ret, err
}

func readByte(file io.Reader) (byte, error) {
	var x byte
	b, err := readNextBytes(file, 1)
	if err != nil {
		return 0, err
	}
	buf := bytes.NewReader(b)
	err = binary.Read(buf, binary.LittleEndian, &x)
	if err != nil {
		return 0, err
	}
	return x, nil
}

func mustReadByte(r io.Reader) byte {
	x, err := readByte(r)
	if err != nil {
		panic(err)
	}
	return x
}

func readString(file io.Reader) (string, error) {
	totalLength := 0
	bt, err := readByte(file)
	if err != nil {
		return "", err
	}
	partialLength := int(bt)
	numBytes := 0
	for partialLength&0x80 > 0 {
		totalLength += int((partialLength & 0x7F) << (7 * numBytes))
		b, err := readByte(file)
		if err != nil {
			return "", err
		}
		partialLength = int(b)
		numBytes++
	}
	totalLength += int(partialLength) << (7 * numBytes)
	b, err := readNextBytes(file, totalLength)
	if err != nil {
		return "", err
	}
	r := string(b)
	if len(r) < totalLength {
		return "", fmt.Errorf("Failed to read complete string")
	}
	return r, nil
}

func mustReadString(r io.Reader) string {
	s, err := readString(r)
	if err != nil {
		panic(err)
	}
	return s
}

func readFloat32(file io.Reader) (float32, error) {
	var r float32
	b, err := readNextBytes(file, 4)
	if err != nil {
		return 0.0, fmt.Errorf("failed to read bytes: %w", err)
	}
	buf := bytes.NewReader(b)
	err = binary.Read(buf, binary.LittleEndian, &r)
	if err != nil {
		return 0.0, fmt.Errorf("failed to convert bytes: %w", err)
	}
	return r, nil

}

func mustReadFloat32(r io.Reader) float32 {
	n, err := readFloat32(r)
	if err != nil {
		panic(err)
	}
	return n
}
