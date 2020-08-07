package beadarray

import (
	"bufio"
	"os"
	"strconv"
	"strings"
)

type CSVBeadPoolManifest struct {
	Names        []string
	LocusEntries map[string]CSVLocusEntry
}

type CSVLocusEntry struct {
	IlmnID          string
	Name            string
	IlmnStrand      string
	SNP             string
	AddressAID      string
	AlleleAProbeSeq string
	AddressBID      string
	AlleleBProbeSeq string
	GenomeBuild     string
	Chr             string
	MapInfo         int
	Ploidy          string
	Species         string
	Source          string
	SourceVersion   string
	SourceStrand    string
	SourceSeq       string
	TopGenomicSeq   string
	BeadSetID       string
	ExpClusters     string
	RefStrand       string
}

func NewCSVBeadPoolManifest(path string) (CSVBeadPoolManifest, error) {
	ret := CSVBeadPoolManifest{}

	ret.LocusEntries = make(map[string]CSVLocusEntry)
	f, err := os.Open(path)
	if err != nil {
		return ret, err
	}
	defer f.Close()
	inData := false
	columnIndexes := make(map[string]int)
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		if line[:6] == "IlmnID" {
			inData = true
			items := strings.Split(line, ",")
			columns := []string{
				"IlmnID",
				"Name",
				"IlmnStrand",
				"SNP",
				"AddressA_ID",
				"AlleleA_ProbeSeq",
				"AddressB_ID",
				"AlleleB_ProbeSeq",
				"GenomeBuild",
				"Chr",
				"MapInfo",
				"Ploidy",
				"Species",
				"Source",
				"SourceVersion",
				"SourceStrand",
				"SourceSeq",
				"TopGenomicSeq",
				"BeadSetID",
				"Exp_Clusters",
				"RefStrand",
			}
			for _, column := range columns {
				columnIndexes[column] = stringSliceIndex(items, column)
			}
			continue
		}
		if line == "[Controls]" {
			inData = false
			break
		}
		if inData {
			items := strings.Split(line, ",")
			name := items[columnIndexes["Name"]]
			ret.Names = append(ret.Names, name)

			mapInfo, err := strconv.Atoi(items[columnIndexes["MapInfo"]])
			if err != nil {
				return ret, err
			}
			// beadsetID
			// expCluster
			ret.LocusEntries[name] = CSVLocusEntry{
				IlmnID:          items[columnIndexes["IlmnID"]],
				Name:            items[columnIndexes["Name"]],
				IlmnStrand:      items[columnIndexes["IlmnStrand"]],
				SNP:             items[columnIndexes["SNP"]],
				AddressAID:      items[columnIndexes["AddressA_ID"]],
				AlleleAProbeSeq: items[columnIndexes["AlleleA_ProbeSeq"]],
				AddressBID:      items[columnIndexes["AddressB_ID"]],
				AlleleBProbeSeq: items[columnIndexes["AlleleB_ProbeSeq"]],
				GenomeBuild:     items[columnIndexes["GenomeBuild"]],
				Chr:             items[columnIndexes["Chr"]],
				MapInfo:         mapInfo,
				Ploidy:          items[columnIndexes["Ploidy"]],
				Species:         items[columnIndexes["Species"]],
				Source:          items[columnIndexes["Source"]],
				SourceVersion:   items[columnIndexes["SourceVersion"]],
				SourceStrand:    items[columnIndexes["SourceStrand"]],
				SourceSeq:       items[columnIndexes["SourceSeq"]],
				TopGenomicSeq:   items[columnIndexes["TopGenomicSeq"]],
				BeadSetID:       items[columnIndexes["BeadSetID"]],
				ExpClusters:     items[columnIndexes["Exp_Clusters"]],
				RefStrand:       items[columnIndexes["RefStrand"]],
			}
		}

	}
	if err := scanner.Err(); err != nil {
		return ret, err
	}
	return ret, nil
}

func stringSliceIndex(list []string, target string) int {
	for i, j := range list {
		if j == target {
			return i
		}
	}
	return -1
}
