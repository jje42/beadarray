package beadarray

import (
	"log"
	"os"
	"reflect"
	"testing"
)

var result float64

// BenchmarkCallRate ...
func BenchmarkCallRate(b *testing.B) {
	fn := "/home/ellisjj/projects/pgx/gsamd-v3/analysis/204379790036/204379790036_R06C01.gtc"
	gtc, err := NewGenotypeCalls(fn)
	if err != nil {
		log.Fatal(err)
	}
	var r float64
	for i := 0; i < b.N; i++ {
		r, _ = gtc.CallRate()
	}
	result = r
}

// func BenchmarkCallRate2(b *testing.B) {
// 	fn := "/home/ellisjj/projects/pgx/gsamd-v3/analysis/204379790036/204379790036_R06C01.gtc"
// 	gtc, err := NewGenotypeCalls(fn)
// 	if err != nil {
// 		log.Fatal(err)
// 	}
// 	var r float64
// 	for i := 0; i < b.N; i++ {
// 		r, _ = gtc.CallRate2()
// 	}
// 	result = r
// }

func TestGenotypeCalls_NormalizedIntensities(t *testing.T) {
	type fields struct {
		file      string
		f         *os.File
		Version   byte
		toc       map[int16]int
		genotypes []byte
	}
	type args struct {
		normalizationLookups []byte
	}
	tests := []struct {
		name    string
		fields  fields
		args    args
		want    []float32
		want1   []float32
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gtc := GenotypeCalls{
				file:      tt.fields.file,
				f:         tt.fields.f,
				Version:   tt.fields.Version,
				toc:       tt.fields.toc,
				genotypes: tt.fields.genotypes,
			}
			got, got1, err := gtc.NormalizedIntensities(tt.args.normalizationLookups)
			if (err != nil) != tt.wantErr {
				t.Errorf("GenotypeCalls.NormalizedIntensities() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("GenotypeCalls.NormalizedIntensities() got = %v, want %v", got, tt.want)
			}
			if !reflect.DeepEqual(got1, tt.want1) {
				t.Errorf("GenotypeCalls.NormalizedIntensities() got1 = %v, want %v", got1, tt.want1)
			}
		})
	}
}
