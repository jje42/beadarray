package beadarray

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"

	"github.com/cheggaaa/pb/v3"
)

// BPM ...
type BPM struct {
	Version          int
	ManifestName     string
	Names            []string
	NumLoci          int
	ControlConfig    string
	NormalizationIDs []byte
	LocusEntries     map[string]LocusEntry
}

// LocusEntry ...
type LocusEntry struct {
	LocusVersion  int
	IlmnID        string
	Name          string
	SNP           string
	Chrom         string
	MapInfo       int
	AddressA      int
	AddressB      int
	AssayType     byte
	RefStrand     string
	GenomeBuild   string
	Source        string
	SourceVersion string
	SourceStrand  string
	Ploidy        string
	Species       string
	IlmnStrand    string
}

//func xmain() {
//path := "/data/manifest/GSAMD/v3/GRCh38/GSAMD-24v3-0-EA_20034606_A2.bpm"
//bpm, err := NewBPM(path)
//if err != nil {
//log.Fatalln(err)
//}
//log.Printf("numLoci = %v", len(bpm.Names))
//log.Printf("%#v", bpm.LocusEntries[bpm.Names[len(bpm.Names)-1]])
//}

// NewBPM ...
func NewBPM(path string) (BPM, error) {
	ret := BPM{}
	f, err := os.Open(path)
	if err != nil {
		return ret, err
	}
	defer f.Close()

	file := bufio.NewReaderSize(f, 4096)

	formatName, err := readNextBytes(file, 3)
	if err != nil {
		return ret, err
	}

	if string(formatName) != "BPM" {
		return ret, fmt.Errorf("file is not BPM format")
	}

	_, err = readNextBytes(file, 1)
	if err != nil {
		return ret, err
	}

	version, err := readInt(file)
	if err != nil {
		return ret, err
	}
	versionFlag := 0x1000

	if (version & versionFlag) == versionFlag {
		version = (version ^ version)
	}
	manifestName, err := readString(file)
	if err != nil {
		return ret, err
	}
	var controlConfig string = ""
	if version > 1 {
		controlConfig, err = readString(file)
		if err != nil {
			return ret, err
		}

	}
	numLoci, err := readInt(file)
	if err != nil {
		return ret, err
	}
	// readNextBytes(file, 4*numLoci)
	for i := 0; i < numLoci; i++ {
		_, err := readInt(file)
		if err != nil {
			return ret, err
		}
	}
	names := make([]string, numLoci)
	for i := 0; i < numLoci; i++ {
		name, err := readString(file)
		if err != nil {
			return ret, err
		}
		names[i] = name
	}
	normalizationIds := make([]byte, numLoci)
	for i := 0; i < numLoci; i++ {
		id, err := readByte(file)
		if err != nil {
			return ret, err
		}
		if id >= 100 {
			return ret, fmt.Errorf("Manifest format error: read invalid normalization ID")
		}
		normalizationIds[i] = id
	}
	locusEntries := make(map[string]LocusEntry)
	bar := pb.StartNew(numLoci)
	for i := 0; i < numLoci; i++ {
		locus, err := NewLocusEntry(file)
		if err != nil {
			return ret, err
		}
		locusEntries[locus.Name] = locus
		bar.Increment()
	}
	bar.Finish()
	return BPM{
		Version:          version,
		ManifestName:     manifestName,
		Names:            names,
		NumLoci:          numLoci,
		ControlConfig:    controlConfig,
		NormalizationIDs: normalizationIds,
		LocusEntries:     locusEntries,
	}, nil
}

// NewLocusEntry ...
func NewLocusEntry(file io.Reader) (LocusEntry, error) {
	ret := LocusEntry{}
	// Read Locus Entry
	locusVersion, err := readInt(file)
	if err != nil {
		return ret, err
	}
	switch locusVersion {
	case 6:
		return ret, fmt.Errorf("can not parse locus entry version 6")
	case 7:
		return ret, fmt.Errorf("can not parse locus entry version 7")
	case 8:
	default:
		return ret, fmt.Errorf("Manifest format error: unknown version for locus entry (%v)", locusVersion)
	}
	ilmnID, err := readString(file)
	if err != nil {
		return ret, err
	}
	name, err := readString(file)
	if err != nil {
		return ret, err
	}

	for j := 0; j < 3; j++ {
		s, err := readString(file)
		if err != nil {
			return ret, err
		}
		if s != "" {
			return ret, fmt.Errorf("a expected empty string, got %s", s)
		}
	}
	// readNextBytes(file, 4)
	// This is a counter from numLoci down...
	_, err = readInt(file)
	if err != nil {
		return ret, err
	}
	// log.Printf("read int = %d\n", i) // 730059

	s, err := readString(file)
	if err != nil {
		return ret, err
	}
	if s != "" {
		return ret, fmt.Errorf("b expected empty string, got %s", s)
	}
	ilmnStrand, err := readString(file)
	if err != nil {
		return ret, err
	}
	snp, err := readString(file)
	if err != nil {
		return ret, err
	}
	chrom, err := readString(file)
	if err != nil {
		return ret, err
	}

	ploidy, err := readString(file)
	if err != nil {
		return ret, err
	}
	species, err := readString(file)
	if err != nil {
		return ret, err
	}

	s, err = readString(file)
	if err != nil {
		return ret, err
	}
	mapInfo, err := strconv.Atoi(s)
	if err != nil {
		return ret, err
	}
	s, err = readString(file)
	if err != nil {
		return ret, err
	}
	if s != "" {
		return ret, fmt.Errorf("c expected empty string, got %s", s)
	}
	sourceStrand, err := readString(file)
	if err != nil {
		return ret, err
	}
	addressA, err := readInt(file)
	if err != nil {
		return ret, err
	}
	addressB, err := readInt(file)
	if err != nil {
		return ret, err
	}

	for i := 0; i < 2; i++ {
		s, err := readString(file)
		if err != nil {
			return ret, err
		}
		if s != "" {
			return ret, fmt.Errorf("d expected empty string, got %s", s)
		}
	}
	genomeBuild, err := readString(file)
	if err != nil {
		return ret, err
	}
	source, err := readString(file)
	if err != nil {
		return ret, err
	}
	sourceVersion, err := readString(file)
	if err != nil {
		return ret, err
	}
	// This appears to be sourceStrand again !?
	_, err = readString(file)
	if err != nil {
		return ret, err
	}
	s, err = readString(file)
	if err != nil {
		return ret, err
	}
	if s != "" {
		return ret, fmt.Errorf("e expected empty string, got %s", s)
	}

	_, err = readNextBytes(file, 3)
	if err != nil {
		return ret, nil
	}
	assayType, err := readByte(file)
	if err != nil {
		return ret, nil
	}

	_, err = readNextBytes(file, 4*4)
	if err != nil {
		return ret, err
	}

	refStrand, err := readString(file)
	if err != nil {
		return ret, err
	}

	return LocusEntry{
		LocusVersion:  locusVersion,
		IlmnID:        ilmnID,
		Name:          name,
		SNP:           snp,
		Chrom:         chrom,
		MapInfo:       mapInfo,
		AddressA:      addressA,
		AddressB:      addressB,
		AssayType:     assayType,
		RefStrand:     refStrand,
		GenomeBuild:   genomeBuild,
		Source:        source,
		SourceVersion: sourceVersion,
		SourceStrand:  sourceStrand,
		Ploidy:        ploidy,
		Species:       species,
		IlmnStrand:    ilmnStrand,
	}, nil
}

// func parseLocusVersion6(f io.Reader) (LocusEntry, error) {
// }

// func parseLocusVersion7(f io.Reader) (LocusEntry, error) {
// }

// func parseLocusVersion8(f io.Reader) (LocusEntry, error) {
// }
