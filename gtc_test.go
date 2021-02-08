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

var nFloats = 730059

func float32ToByte(f float32) []byte {
	var buf [4]byte
	binary.LittleEndian.PutUint32(buf[:], math.Float32bits(f))
	return buf[:]
}

func readFloat32Slice1(r io.Reader, n int) []float32 {
	xs := make([]float32, nFloats)
	err := binary.Read(r, binary.LittleEndian, &xs)
	if err != nil {
		panic(err)
	}
	return xs
}

func readFloat32Slice2(r io.Reader, n int) []float32 {
	bs := make([]byte, int(n)*4)
	if _, err := io.ReadAtLeast(r, bs, len(bs)); err != nil {
		panic(err)
	}
	xs := make([]float32, n)
	for i := 0; i < int(n); i++ {
		start := i * 4
		end := start + 4
		x := math.Float32frombits(binary.LittleEndian.Uint32(bs[start:end]))
		xs[i] = x
	}
	return xs
}

func readFloat32Slice3(r io.Reader, n int) []float32 {
	buf := make([]byte, int(n)*4)
	if _, err := io.ReadAtLeast(r, buf, len(buf)); err != nil {
		log.Fatal(err)
	}
	header := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
	header.Len /= 4
	header.Cap /= 4
	return *(*[]float32)(unsafe.Pointer(&header))
}
func setup() ([]byte, []float32) {
	bs := []byte{}
	expected := []float32{}
	for i := 0; i < nFloats; i++ {
		bs = append(bs, float32ToByte(float32(i))...)
		expected = append(expected, float32(i))
	}
	return bs, expected
}

func BenchmarkTestIt1(b *testing.B) {
	bs, _ := setup()
	for i := 0; i < b.N; i++ {
		_ = readFloat32Slice1(bytes.NewReader(bs), nFloats)
	}
}

func BenchmarkTestIt2(b *testing.B) {
	bs, _ := setup()
	for i := 0; i < b.N; i++ {
		_ = readFloat32Slice2(bytes.NewReader(bs), nFloats)
	}
}
func BenchmarkTestIt3(b *testing.B) {
	bs, _ := setup()
	for i := 0; i < b.N; i++ {
		_ = readFloat32Slice3(bytes.NewReader(bs), nFloats)
	}
}

func Test_readFloat32Slice1(t *testing.T) {
	bs, expected := setup()
	r := bytes.NewReader(bs)
	type args struct {
		r io.Reader
		n int
	}
	tests := []struct {
		name string
		args args
		want []float32
	}{
		{"t1", args{r, nFloats}, expected},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := readFloat32Slice1(tt.args.r, tt.args.n); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("readFloat32Slice1() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_readFloat32Slice2(t *testing.T) {
	bs, expected := setup()
	r := bytes.NewReader(bs)
	type args struct {
		r io.Reader
		n int
	}
	tests := []struct {
		name string
		args args
		want []float32
	}{
		{"t1", args{r, nFloats}, expected},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := readFloat32Slice2(tt.args.r, tt.args.n); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("readFloat32Slice2() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_readFloat32Slice3(t *testing.T) {
	bs, expected := setup()
	r := bytes.NewReader(bs)
	type args struct {
		r io.Reader
		n int
	}
	tests := []struct {
		name string
		args args
		want []float32
	}{
		{"t1", args{r, nFloats}, expected},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := readFloat32Slice3(tt.args.r, tt.args.n); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("readFloat32Slice3() = %v, want %v", got, tt.want)
			}
		})
	}
}
