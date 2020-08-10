package beadarray

import (
	"encoding/binary"
	"fmt"
	"io"
)

// EGT ...
type EGT struct {
	GencallVersion       string
	ClusterVersion       string
	CallVersion          string
	NormalizationVersion string
	DateCreated          string
	ManifestName         string
	Name2ClusterRecord   map[string]ClusterRecord
}

// ClusterRecord  ...
type ClusterRecord struct {
	AAClusterStats     ClusterStats
	ABClusterStats     ClusterStats
	BBClusterStats     ClusterStats
	IntensityThreshold float32
	ClusterScore       ClusterScore
	Address            int
}

// ClusterScore ...
type ClusterScore struct {
	ClusterSeparation float32
	TotalScore        float32
	OriginalScore     float32
	Edited            bool
}

// ClusterStats ...
type ClusterStats struct {
	ThetaMean float32
	ThetaDev  float32
	RMean     float32
	RDev      float32
	N         int
}

// NewEGT ...
func NewEGT(r io.Reader) (egt *EGT, err error) {
	egt = &EGT{}
	defer func() {
		if r := recover(); r != nil {
			err = r.(error)
		}
	}()
	version := mustReadInt(r)
	if version != 3 {
		return nil, fmt.Errorf("Cluster file version %d not supported", version)
	}
	gencallVersion := mustReadString(r)
	clusterVersion := mustReadString(r)
	callVersion := mustReadString(r)
	normalizationVersion := mustReadString(r)
	dateCreated := mustReadString(r)
	isWgt := mustReadByte(r)
	if isWgt == 0 {
		return nil, fmt.Errorf("Only WGT cluster file version supported")
	}
	manifestName := mustReadString(r)
	dataBlockVersion := mustReadInt(r)
	if dataBlockVersion != 8 && dataBlockVersion != 9 {
		return nil, fmt.Errorf("Data block version in cluster file %d not supported", dataBlockVersion)
	}

	//opa
	mustReadString(r) // file name?

	numRecords := mustReadInt(r)
	clusterRecords := make([]ClusterRecord, numRecords)
	// won't work as there is extra data
	// err = binary.Read(r, binary.LittleEndian, &clusterRecords)
	// if err != nil {
	// return nil, err
	// }
	for i := 0; i < numRecords; i++ {
		record, err := readClusterRecord(r, dataBlockVersion)
		if err != nil {
			return nil, err
		}
		clusterRecords[i] = record
	}

	clusterScores := make([]ClusterScore, numRecords)
	for i := 0; i < numRecords; i++ {
		clusterScores[i] = readClusterScore(r)
	}
	// genotypes
	for i := 0; i < numRecords; i++ {
		mustReadString(r)
	}
	lociNames := make([]string, numRecords)
	for i := 0; i < numRecords; i++ {
		lociNames[i] = mustReadString(r)
	}

	addresses := make([]int32, numRecords)
	// for i := 0; i < numRecords; i++ {
	// 	addresses[i] = mustReadInt(r)
	// }
	// Is this "better"?
	if err := binary.Read(r, binary.LittleEndian, &addresses); err != nil {
		return egt, err
	}

	// cluster counts. In BeadArrayFiles these counts are only used to
	// assert they are equal to the Ns in the cluster records.
	for i := 0; i < numRecords; i++ {
		mustReadInt(r) // AA
		mustReadInt(r) // AB
		mustReadInt(r) // BB
	}

	//  Add address and cluster_score to each record.
	name2clusterRecord := make(map[string]ClusterRecord)
	for i := 0; i < numRecords; i++ {
		record := clusterRecords[i]
		record.Address = int(addresses[i])
		record.ClusterScore = clusterScores[i]
		name2clusterRecord[lociNames[i]] = record
	}

	egt.GencallVersion = gencallVersion
	egt.ClusterVersion = clusterVersion
	egt.CallVersion = callVersion
	egt.NormalizationVersion = normalizationVersion
	egt.DateCreated = dateCreated
	egt.ManifestName = manifestName
	egt.Name2ClusterRecord = name2clusterRecord
	return egt, nil
}

func readClusterRecord(r io.Reader, version int) (ClusterRecord, error) {
	record := ClusterRecord{}
	if version != 9 {
		return record, fmt.Errorf("unsupported cluster record version %d", version)
	}
	aaN := mustReadInt(r)
	abN := mustReadInt(r)
	bbN := mustReadInt(r)
	ys := make([]float32, 13)
	if err := binary.Read(r, binary.LittleEndian, &ys); err != nil {
		return record, err
	}
	// aaRDev := mustReadFloat32(r)
	// abRDev := mustReadFloat32(r)
	// bbRDev := mustReadFloat32(r)
	// aaRMean := mustReadFloat32(r)
	// abRMean := mustReadFloat32(r)
	// bbRMean := mustReadFloat32(r)
	// aaThetaDev := mustReadFloat32(r)
	// abThetaDev := mustReadFloat32(r)
	// bbThetaDev := mustReadFloat32(r)
	// aaThetaMean := mustReadFloat32(r)
	// abThetaMean := mustReadFloat32(r)
	// bbThetaMean := mustReadFloat32(r)
	// intensityThreshold := mustReadFloat32(r)

	// read through unused fields
	for i := 0; i < 14; i++ {
		mustReadFloat32(r)
	}
	return ClusterRecord{
		AAClusterStats: ClusterStats{
			ThetaMean: ys[9], // aaThetaMean,
			ThetaDev:  ys[6], // aaThetaDev,
			RMean:     ys[3], // aaRMean,
			RDev:      ys[0], // aaRDev,
			N:         aaN,   // aaN,
		},
		ABClusterStats: ClusterStats{
			ThetaMean: ys[10], // abThetaMean,
			ThetaDev:  ys[7],  // abThetaDev,
			RMean:     ys[4],  // abRMean,
			RDev:      ys[1],  // abRDev,
			N:         abN,    // abN,
		},
		BBClusterStats: ClusterStats{
			ThetaMean: ys[11], // bbThetaMean,
			ThetaDev:  ys[8],  // bbThetaDev,
			RMean:     ys[5],  // bbRMean,
			RDev:      ys[2],  // bbRDev,
			N:         bbN,    // bbN,
		},
		IntensityThreshold: ys[12], //intensityThreshold,
	}, nil
}

func readClusterScore(r io.Reader) ClusterScore {
	return ClusterScore{
		ClusterSeparation: mustReadFloat32(r),
		TotalScore:        mustReadFloat32(r),
		OriginalScore:     mustReadFloat32(r),
		Edited:            mustReadByte(r) != 0,
	}
}
