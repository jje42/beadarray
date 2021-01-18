package beadarray

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
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
func NewBPM(path string) (ret BPM, err error) {
	ret = BPM{}
	defer func() {
		if r := recover(); r != nil {
			err = r.(error)
		}
	}()

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

	version := mustReadInt(file)
	versionFlag := 0x1000

	if (version & versionFlag) == versionFlag {
		version = (version ^ version)
	}
	manifestName := mustReadString(file)
	var controlConfig string = ""
	if version > 1 {
		controlConfig = mustReadString(file)
	}
	numLoci := mustReadInt(file)
	// readNextBytes(file, 4*numLoci)
	for i := 0; i < numLoci; i++ {
		_, err := readInt(file)
		if err != nil {
			return ret, err
		}
	}
	names := make([]string, numLoci)
	for i := 0; i < numLoci; i++ {
		names[i] = mustReadString(file)
	}
	normalizationIds := make([]byte, numLoci)
	for i := 0; i < numLoci; i++ {
		id := mustReadByte(file)
		if id >= 100 {
			return ret, fmt.Errorf("Manifest format error: read invalid normalization ID")
		}
		normalizationIds[i] = id
	}
	locusEntries := make(map[string]LocusEntry)
	for i := 0; i < numLoci; i++ {
		locus, err := NewLocusEntry(file)
		if err != nil {
			return ret, err
		}
		locusEntries[locus.Name] = locus
	}
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
	locusVersion := mustReadInt(file)
	switch locusVersion {
	case 6:
		return ret, fmt.Errorf("can not parse locus entry version 6")
	case 7:
		return ret, fmt.Errorf("can not parse locus entry version 7")
	case 8:
	default:
		return ret, fmt.Errorf("Manifest format error: unknown version for locus entry (%v)", locusVersion)
	}
	ilmnID := mustReadString(file)
	name := mustReadString(file)

	for j := 0; j < 3; j++ {
		s := mustReadString(file)
		if s != "" {
			return ret, fmt.Errorf("a expected empty string, got %s", s)
		}
	}
	// This is a counter from numLoci down...
	_ = mustReadInt(file)

	s := mustReadString(file)
	if s != "" {
		return ret, fmt.Errorf("b expected empty string, got %s", s)
	}
	ilmnStrand := mustReadString(file)
	snp := mustReadString(file)
	chrom := mustReadString(file)
	ploidy := mustReadString(file)
	species := mustReadString(file)
	s = mustReadString(file)
	mapInfo, err := strconv.Atoi(s)
	if err != nil {
		return ret, err
	}
	s = mustReadString(file)
	if s != "" {
		return ret, fmt.Errorf("c expected empty string, got %s", s)
	}
	sourceStrand := mustReadString(file)
	addressA := mustReadInt(file)
	addressB := mustReadInt(file)

	for i := 0; i < 2; i++ {
		s := mustReadString(file)
		if s != "" {
			return ret, fmt.Errorf("d expected empty string, got %s", s)
		}
	}
	genomeBuild := mustReadString(file)
	source := mustReadString(file)
	sourceVersion := mustReadString(file)
	// This appears to be sourceStrand again !?
	_ = mustReadString(file)
	s = mustReadString(file)
	if s != "" {
		return ret, fmt.Errorf("e expected empty string, got %s", s)
	}

	_, err = readNextBytes(file, 3)
	if err != nil {
		return ret, nil
	}
	assayType := mustReadByte(file)

	_, err = readNextBytes(file, 4*4)
	if err != nil {
		return ret, err
	}
	refStrand := mustReadString(file)

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
