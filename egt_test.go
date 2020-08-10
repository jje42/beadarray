package beadarray

import (
	"bytes"
	"encoding/binary"
	"io"
	"log"
	"math"
	"reflect"
	"testing"
	"unsafe"
)

var n = 730059
var bb bytes.Buffer

func init() {
	bb = bytes.Buffer{}
	for i := 0; i < n; i++ {
		buf := make([]byte, 4)
		binary.LittleEndian.PutUint32(buf, uint32(i))
		_, err := bb.Write(buf)
		if err != nil {
			log.Fatal(err)
		}
	}
}

func BenchmarkReadSlice(b *testing.B) {
	for i := 0; i < b.N; i++ {
		r := bytes.NewReader(bb.Bytes())
		x := make([]int32, n)
		binary.Read(r, binary.LittleEndian, &x)
	}
}

func BenchmarkReadMustReadInt(b *testing.B) {
	for i := 0; i < b.N; i++ {
		r := bytes.NewReader(bb.Bytes())
		xs := make([]int32, n)
		for j := 0; j < n; j++ {
			xs[j] = int32(mustReadInt(r))
		}
	}
}

func BenchmarkReadUnsafe(b *testing.B) {
	for i := 0; i < b.N; i++ {
		r := bytes.NewReader(bb.Bytes())
		buf := make([]byte, n*4)
		if _, err := io.ReadAtLeast(r, buf, len(buf)); err != nil {
			panic(err)
		}
		header := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
		header.Len /= 4
		header.Cap /= 4
		_ = *(*[]int32)(unsafe.Pointer(&header))
	}
}

func BenchmarkReadByteOrder(b *testing.B) {
	for i := 0; i < b.N; i++ {
		xs := make([]int32, n)
		for j := 0; j < n; j++ {
			start := j * 4
			end := j + 4
			x := binary.LittleEndian.Uint32(bb.Bytes()[start:end])
			xs[j] = int32(x)
		}
	}
}

func BenchmarkReadFloat32ByteOrder(b *testing.B) {
	for i := 0; i < b.N; i++ {
		xs := make([]float32, n)
		for j := 0; j < n; j++ {
			start := j * 4
			end := start + 4
			x := binary.LittleEndian.Uint32(bb.Bytes()[start:end])
			xs[j] = math.Float32frombits(x)
		}
	}
}

func BenchmarkReadFloat32Unsafe(b *testing.B) {
	for i := 0; i < b.N; i++ {
		r := bytes.NewReader(bb.Bytes())
		buf := make([]byte, n*4)
		if _, err := io.ReadAtLeast(r, buf, len(buf)); err != nil {
			panic(err)
		}
		header := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
		header.Len /= 4
		header.Cap /= 4
		_ = *(*[]int32)(unsafe.Pointer(&header))
	}
}

func TestReadTest1(t *testing.T) {
	r := bytes.NewReader(bb.Bytes())
	xs := make([]int32, n)
	binary.Read(r, binary.LittleEndian, &xs)
	if len(xs) != n {
		t.Errorf("returned array length not equal to %v, got %v", n, len(xs))
	}
}

func TestReadTest2(t *testing.T) {
	r := bytes.NewReader(bb.Bytes())
	xs := make([]int32, n)
	for j := 0; j < n; j++ {
		xs[j] = int32(mustReadInt(r))
	}

	if len(xs) != n {
		t.Errorf("returned array length not equal to %v, got %v", n, len(xs))
	}
}

func TestReadTest3(t *testing.T) {
	r := bytes.NewReader(bb.Bytes())
	buf := make([]byte, n*4)
	if _, err := io.ReadAtLeast(r, buf, len(buf)); err != nil {
		panic(err)
	}
	header := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
	header.Len /= 4
	header.Cap /= 4
	xs := *(*[]int32)(unsafe.Pointer(&header))
	if len(xs) != n {
		t.Errorf("returned array length not equal to %v, got %v", n, len(xs))
	}
}
