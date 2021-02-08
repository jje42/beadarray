package beadarray

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"net/http"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"unsafe"
)

// GTC ...
type GTC struct {
	file      string
	f         *os.File // bufio.Reader has no Seek method
	Version   byte
	toc       map[int16]int
	genotypes []byte
}

const (
	idNumSnps                 = 1
	idPloidy                  = 2
	idPloidyType              = 3
	idSampleName              = 10
	idSamplePlate             = 11
	idSampleWell              = 12
	idClusterFile             = 100
	idSnpManifest             = 101
	idImagingDate             = 200
	idAutocallDate            = 201
	idAutocallVersion         = 300
	idNormalizationTransforms = 400
	idControlsX               = 500
	idControlsY               = 501
	idRawX                    = 1000
	idRawY                    = 1001
	idGenotypes               = 1002
	idBaseCalls               = 1003
	idGenotypeScores          = 1004
	idScannerData             = 1005
	idCallRate                = 1006
	idGender                  = 1007
	idLogrDev                 = 1008
	idGc10                    = 1009
	idGc50                    = 1011
	idBAlleleFreqs            = 1012
	idLogrRatios              = 1013
	idPercentilesX            = 1014
	idPercentilesY            = 1015
	idSlideIdentifier         = 1016
)

var code2genotype = []string{
	"NC",
	"AA",
	"AB",
	"BB",
	"NULL",
	"A",
	"B",
	"AAA",
	"AAB",
	"ABB",
	"BBB",
	"AAAA",
	"AAAB",
	"AABB",
	"ABBB",
	"BBBB",
	"AAAAA",
	"AAAAB",
	"AAABB",
	"AABBB",
	"ABBBB",
	"BBBBB",
	"AAAAAA",
	"AAAAAB",
	"AAAABB",
	"AAABBB",
	"AABBBB",
	"ABBBBB",
	"BBBBBB",
	"AAAAAAA",
	"AAAAAAB",
	"AAAAABB",
	"AAAABBB",
	"AAABBBB",
	"AABBBBB",
	"ABBBBBB",
	"BBBBBBB",
	"AAAAAAAA",
	"AAAAAAAB",
	"AAAAAABB",
	"AAAAABBB",
	"AAAABBBB",
	"AAABBBBB",
	"AABBBBBB",
	"ABBBBBBB",
	"BBBBBBBB",
}

// Code2Genotype ...
var Code2Genotype = code2genotype

// IsGTCFile returns true if file appears to be a valid GTC file and false
// otherwise.
func IsGTCFile(file string) (bool, error) {
	// Content type detection is used to avoid falsely returning true when a
	// text file is passed that happens to have "gtc" in its first three
	// bytes. For example, this can happen in FOFNs when the GTCs are in a
	// directory named "gtc".
	r, err := os.Open(file)
	if err != nil {
		return false, fmt.Errorf("failed to open file: %w", err)
	}
	defer r.Close()
	bs := make([]byte, 512)
	n, err := r.Read(bs)
	if err != nil {
		return false, fmt.Errorf("unable to read file: %w", err)
	}
	contentType := http.DetectContentType(bs[:n])
	return contentType == "application/octet-stream" && n > 3 && string(bs[:3]) == "gtc", nil
}

// NewGTC ...
func NewGTC(file string) (GTC, error) {
	f, err := os.Open(file)
	if err != nil {
		return GTC{}, err
	}

	identifier, err := readNextBytes(f, 3)
	if err != nil {
		return GTC{}, err
	}
	if string(identifier) != "gtc" {
		return GTC{}, fmt.Errorf("GTC format error: bad format identifier")
	}
	version, err := readByte(f)
	if err != nil {
		return GTC{}, err
	}
	if !isSupportedVersion(version) {
		return GTC{}, fmt.Errorf("Unsupported GTF file version (%v)", version)
	}

	n, err := readInt(f)
	if err != nil {
		return GTC{}, err
	}

	toc := make(map[int16]int)

	for i := 0; i < n; i++ {
		id, err := readInt16(f)
		if err != nil {
			return GTC{}, err
		}
		offset, err := readInt(f)
		if err != nil {
			return GTC{}, err
		}
		toc[id] = offset
	}

	return GTC{file: file, f: f, Version: version, toc: toc}, nil
}

// Filename ...
func (g GTC) Filename() string {
	return g.file
}

// CallRate ...
func (g GTC) CallRate() (float32, error) {
	return g.genericFloat(idCallRate)
}

// LogRDev ...
func (g GTC) LogRDev() (float32, error) {
	return g.genericFloat(idLogrDev)
}

// GC10 returns the GC10 (GenCall score - 10th percentile).
func (g GTC) GC10() (float32, error) {
	return g.genericFloat(idGc10)
}

// GC50 returns the GC50 (GenCall score - 50th percentile).
func (g GTC) GC50() (float32, error) {
	return g.genericFloat(idGc50)
}

// NumCalls returns the number of calls.
func (g GTC) NumCalls() (int, error) {
	return g.genericInt(g.toc[idGc50] + 4)
}

// NumNoCalls ...
func (g GTC) NumNoCalls() (int, error) {
	return g.genericInt(g.toc[idGc50] + 8)
}

// NumIntensityOnly returns the number of intensity only SNPs
func (g GTC) NumIntensityOnly() (int, error) {
	return g.genericInt(g.toc[idGc50] + 12)
}

// SampleName ...
func (g GTC) SampleName() (string, error) {
	name, err := g.genericString(g.toc[idSampleName])
	if err != nil {
		return "", err
	}
	if name == "" {
		return strings.TrimSuffix(filepath.Base(g.file), ".gtc"), nil
	}
	return name, nil
}

// ClusterFile ...
func (g GTC) ClusterFile() (string, error) {
	return g.genericString(g.toc[idClusterFile])
}

// SlideIdentifier ...
func (g GTC) SlideIdentifier() (string, error) {
	return g.genericString(g.toc[idSlideIdentifier])
}

// SamplePlate ...
func (g GTC) SamplePlate() (string, error) {
	return g.genericString(g.toc[idSamplePlate])
}

// SampleWell ...
func (g GTC) SampleWell() (string, error) {
	return g.genericString(g.toc[idSampleWell])
}

// SnpManifest ...
func (g GTC) SnpManifest() (string, error) {
	return g.genericString(g.toc[idSnpManifest])
}

// ImagingDate ...
func (g GTC) ImagingDate() (string, error) {
	return g.genericString(g.toc[idImagingDate])
}

// AutocallDate ...
func (g GTC) AutocallDate() (string, error) {
	return g.genericString(g.toc[idAutocallDate])
}

// AutocallVersion ...
func (g GTC) AutocallVersion() (string, error) {
	return g.genericString(g.toc[idAutocallVersion])
}

// Gender ...
func (g GTC) Gender() (string, error) {
	pos := g.toc[idGender]
	_, err := g.f.Seek(int64(pos), io.SeekStart)
	if err != nil {
		return "", err
	}
	r, err := readByte(g.f)
	if err != nil {
		return "", err
	}
	return string(r), nil
}

// Genotypes ...
func (g *GTC) Genotypes() ([]byte, error) {
	if g.genotypes != nil {
		// Just return cached result if available
		return g.genotypes, nil
	}

	// f, err := os.Open(gtc.file)
	// if err != nil {
	// 	return nil, err
	// }
	// defer f.Close()
	pos := g.toc[idGenotypes]
	_, err := g.f.Seek(int64(pos), io.SeekStart)
	if err != nil {
		return nil, err
	}
	b := bufio.NewReader(g.f)
	numEntries, err := readInt(b)
	if err != nil {
		return nil, err
	}
	r := make([]byte, numEntries)
	for i := 0; i < numEntries; i++ {
		x, err := readByte(b)
		if err != nil {
			return nil, err
		}
		r[i] = x
	}
	// Cache result
	g.genotypes = r
	return r, nil
}

// PloidyType ...
func (g GTC) PloidyType() int {
	return g.toc[idPloidyType]
}

// BaseCalls ...
func (g GTC) BaseCalls() ([]string, error) {
	ploidyType := g.PloidyType()
	genotypes, err := g.Genotypes()
	if err != nil {
		return nil, err
	}
	// f, err := os.Open(gtc.file)
	// if err != nil {
	// 	return nil, err
	// }
	// defer f.Close()
	pos := int64(g.toc[idBaseCalls])
	_, err = g.f.Seek(pos, io.SeekStart)
	if err != nil {
		return nil, err
	}
	b := bufio.NewReader(g.f)
	numEntries, err := readInt(b)
	if err != nil {
		return nil, err
	}
	r := make([]string, numEntries)
	for i := 0; i < numEntries; i++ {
		if ploidyType == 1 {
			bytes, err := readNextBytes(b, 2)
			if err != nil {
				return nil, err
			}
			r[i] = string(bytes)
		} else {
			bytes, err := readNextBytes(b, 2)
			if err != nil {
				return nil, err
			}
			byteString := string(bytes)
			abGenotype := code2genotype[genotypes[i]]
			if abGenotype == "NC" || abGenotype == "NULL" {
				r[i] = "-"
			} else {
				var topGenotype []byte
				for j := 0; j < len(abGenotype); j++ {
					if abGenotype[j] == 'A' {
						topGenotype = append(topGenotype, byteString[0])
					} else {
						topGenotype = append(topGenotype, byteString[1])
					}
					r[i] = string(topGenotype)
				}
			}

		}
	}
	return r, nil
}

// Close ...
func (g GTC) Close() {
	g.f.Close()
}

// BAlleleFreqs ...
func (g GTC) BAlleleFreqs() ([]float32, error) {
	if g.Version < 4 {
		return nil, fmt.Errorf("B allele frequencies unavailable in GTC File version %v", g.Version)
	}
	return g.genericFloat32Slice(idBAlleleFreqs)
}

// LogRRatios ...
func (g GTC) LogRRatios() ([]float32, error) {
	if g.Version < 4 {
		return nil, fmt.Errorf("LogR ratios unavailable in GTC File version %v", g.Version)
	}
	return g.genericFloat32Slice(idLogrRatios)
}

// NormalizationTransform ...
type NormalizationTransform struct {
	Version int
	OffsetX float32
	OffsetY float32
	ScaleX  float32
	ScaleY  float32
	Shear   float32
	Theta   float32
}

func readNormalizationTransform(r io.Reader) (NormalizationTransform, error) {
	ret := NormalizationTransform{}
	version, err := readInt(r)
	if err != nil {
		return ret, err
	}
	offsetX, err := readFloat32(r)
	if err != nil {
		return ret, err
	}
	offsetY, err := readFloat32(r)
	if err != nil {
		return ret, err
	}
	scaleX, err := readFloat32(r)
	if err != nil {
		return ret, err
	}
	scaleY, err := readFloat32(r)
	if err != nil {
		return ret, err
	}
	shear, err := readFloat32(r)
	if err != nil {
		return ret, err
	}
	theta, err := readFloat32(r)
	if err != nil {
		return ret, err
	}
	return NormalizationTransform{
		Version: version,
		OffsetX: offsetX,
		OffsetY: offsetY,
		ScaleX:  scaleX,
		ScaleY:  scaleY,
		Shear:   shear,
		Theta:   theta,
	}, nil
}

// RectToPolar ...
func (nt NormalizationTransform) RectToPolar(x, y float32) (float64, float64) {
	if x == 0 && y == 0 {
		return math.NaN(), math.NaN()
	}
	return float64(x + y), math.Atan2(float64(y), float64(x)) * 2.0 / math.Pi
}

// NormalizeIntensities ...
func (nt NormalizationTransform) NormalizeIntensities(x, y float32, threshold bool) (float32, float32) {
	if x == 0 && y == 0 {
		return float32(math.NaN()), float32(math.NaN())
	}
	theta := float64(nt.Theta)

	tempx := float64(x - nt.OffsetX)
	tempy := float64(y - nt.OffsetY)

	tempx2 := math.Cos(theta)*tempx + math.Sin(theta)*tempy
	tempy2 := -math.Sin(theta)*tempx + math.Cos(theta)*tempy

	tempx3 := tempx2 - float64(nt.Shear)*tempy2
	tempy3 := tempy2

	xn := tempx3 / float64(nt.ScaleX)
	yn := tempy3 / float64(nt.ScaleY)

	if threshold {
		if 0 > xn {
			xn = 0
		}
		if 0 > yn {
			yn = 0
		}
	}
	return float32(xn), float32(yn)
}

// NormalizationTransforms ...
func (g GTC) NormalizationTransforms() ([]NormalizationTransform, error) {
	// The components of a NormalizationTransform do not sum to 52 bytes,
	// but they are in 52 byte blocks. Must read in 52 bytes then extract
	// from that. I don't know what, if anything, is in the remaining bytes.
	pos := g.toc[idNormalizationTransforms]
	_, err := g.f.Seek(int64(pos), io.SeekStart)
	if err != nil {
		return nil, fmt.Errorf("GTC.NormalizationTransforms failed: %w", err)
	}
	b := bufio.NewReader(g.f)
	numEntries, _ := readInt(b)
	r := make([]NormalizationTransform, numEntries)
	for i := 0; i < numEntries; i++ {
		buf, err := readNextBytes(b, 52)
		if err != nil {
			return nil, fmt.Errorf("GTC.NormalizationTransforms failed: %w", err)
		}
		nt, err := readNormalizationTransform(bytes.NewReader(buf))
		if err != nil {
			return nil, fmt.Errorf("GTC.NormalizationTransforms failed: %w", err)
		}
		r[i] = nt
	}
	return r, nil
}

// NormalizedIntensities ...
func (g GTC) NormalizedIntensities(normalizationLookups []byte) ([]float32, []float32, error) {
	normalizationTransforms, err := g.NormalizationTransforms()
	if err != nil {
		return nil, nil, err
	}
	rawX, err := g.RawXIntensities()
	rawY, err := g.RawYIntensities()

	xs := make([]float32, len(rawX))
	ys := make([]float32, len(rawY))
	for i := 0; i < len(rawX); i++ {
		x := rawX[i]
		y := rawY[i]
		lookup := normalizationLookups[i]
		// nt := normalizationTransforms[lookup]
		nx, ny := normalizationTransforms[lookup].NormalizeIntensities(float32(x), float32(y), true)
		xs[i] = nx
		ys[i] = ny
		// log.Printf("rawX = %v, rawY = %v, nx = %v, ny = %v, nt = %#v\n", x, y, nx, ny, nt)
	}
	return xs, ys, nil
}

// GenotypeScores returns the genotype scores.
func (g GTC) GenotypeScores() ([]float32, error) {
	return g.genericFloat32Slice(idGenotypeScores)
}

// ControlXIntensities returns the x intensities of control bead types.
func (g GTC) ControlXIntensities() ([]uint16, error) {
	return g.genericUint16Slice(idControlsX)
}

// ControlYIntensities returns the y intensities of control bead types.
func (g GTC) ControlYIntensities() ([]uint16, error) {
	return g.genericUint16Slice(idControlsY)
}

// RawXIntensities returns the raw x intensities of assay bead types.
func (g GTC) RawXIntensities() ([]int16, error) {
	return g.genericInt16Slice(idRawX)
}

// RawYIntensities returns the raw y intensities of assay bead types.
func (g GTC) RawYIntensities() ([]int16, error) {
	return g.genericInt16Slice(idRawY)
}

// ScannerData ...
type ScannerData struct {
	Name     string
	PmtGreen int
	PmtRed   int
	Version  string
	User     string
}

func readScannerData(r io.Reader) ScannerData {
	ret := ScannerData{}
	ret.Name = mustReadString(r)
	ret.PmtGreen = mustReadInt(r)
	ret.PmtRed = mustReadInt(r)
	ret.Version = mustReadString(r)
	ret.User = mustReadString(r)
	return ret
}

// ScannerData returns information about scanner
func (g GTC) ScannerData() ScannerData {
	pos := g.toc[idScannerData]
	_, _ = g.f.Seek(int64(pos), io.SeekStart)
	b := bufio.NewReader(g.f)
	return readScannerData(b)
}

// PercentilesX returns a slice of length three representing 5th, 50th and 95th percentile for x intensity.
func (g GTC) PercentilesX() ([]uint16, error) {
	return g.genericUint16Slice(idPercentilesX)
}

// PercentileY returns a slice of length three representing 5th, 50th and 95th percentile for y intensity.
func (g GTC) PercentileY() ([]uint16, error) {
	return g.genericUint16Slice(idPercentilesY)
}

func (g GTC) gotoPosition(tocEntry int16) (io.Reader, error) {
	pos := int64(g.toc[tocEntry])
	_, err := g.f.Seek(pos, io.SeekStart)
	if err != nil {
		return nil, err
	}
	b := bufio.NewReader(g.f)
	return b, nil
}

func (g GTC) genericString(tocEntry int) (string, error) {
	_, err := g.f.Seek(int64(tocEntry), io.SeekStart)
	if err != nil {
		return "", err
	}
	return readString(g.f)
}

func (g GTC) genericInt(pos int) (int, error) {
	_, err := g.f.Seek(int64(pos), io.SeekStart)
	if err != nil {
		return 0, err
	}
	return readInt(g.f)
}

func (g GTC) genericFloat(tocEntry int16) (float32, error) {
	b, err := g.gotoPosition(tocEntry)
	if err != nil {
		return 0.0, err
	}
	var r float32
	err = binary.Read(b, binary.LittleEndian, &r)
	return r, err
}

func (g *GTC) genericInt16Slice(tocEntry int16) ([]int16, error) {
	b, err := g.gotoPosition(tocEntry)
	if err != nil {
		return nil, err
	}
	var numEntries int32
	if err := binary.Read(b, binary.LittleEndian, &numEntries); err != nil {
		return nil, err
	}
	buf := make([]byte, int(numEntries)*2)
	if _, err := io.ReadAtLeast(b, buf, len(buf)); err != nil {
		return nil, err
	}
	header := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
	header.Len /= 2
	header.Cap /= 2
	xs := *(*[]int16)(unsafe.Pointer(&header))
	return xs, nil
}

func (g *GTC) genericUint16Slice(tocEntry int16) ([]uint16, error) {
	b, err := g.gotoPosition(tocEntry)
	if err != nil {
		return nil, err
	}
	var numEntries int32
	if err := binary.Read(b, binary.LittleEndian, &numEntries); err != nil {
		return nil, err
	}
	buf := make([]byte, int(numEntries)*2)
	if _, err := io.ReadAtLeast(b, buf, len(buf)); err != nil {
		return nil, err
	}
	header := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
	header.Len /= 2
	header.Cap /= 2
	xs := *(*[]uint16)(unsafe.Pointer(&header))
	return xs, nil
}

func (g *GTC) genericFloat32Slice(tocEntry int16) ([]float32, error) {
	b, err := g.gotoPosition(tocEntry)
	if err != nil {
		return nil, err
	}
	var numEntries int32
	if err := binary.Read(b, binary.LittleEndian, &numEntries); err != nil {
		return nil, err
	}
	buf := make([]byte, int(numEntries)*4)
	if _, err := io.ReadAtLeast(b, buf, len(buf)); err != nil {
		return nil, err
	}
	header := *(*reflect.SliceHeader)(unsafe.Pointer(&buf))
	header.Len /= 4
	header.Cap /= 4
	xs := *(*[]float32)(unsafe.Pointer(&header))

	// math.Float32frombits uses unsafe anyway so why would this approach be
	// "better" or "safer"?
	// xs := make([]float32, numEntries)
	// for i := 0; i < int(numEntries); i++ {
	// 	start := i * 4
	// 	end := start + 4
	// 	x := math.Float32frombits(binary.LittleEndian.Uint32(buf[start:end]))
	// 	xs[i] = x
	// }

	return xs, nil
}

func isSupportedVersion(version byte) bool {
	r := false
	supportedVersions := []byte{3, 4, 5}
	for _, value := range supportedVersions {
		if value == version {
			r = true
		}
	}
	return r
}
