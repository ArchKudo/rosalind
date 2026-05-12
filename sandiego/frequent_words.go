package main

import (
	"fmt"
)

func main() {

	hay := "CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA"
	k := 3
	// var freq map[string]int = {}
	freq := make(map[string]int, len(hay))

	for i := 0; i < len(hay)-k+1; i++ {
		kmer := hay[i : i+k]

		if _, ok := freq[kmer]; ok {
			freq[kmer] += 1
		} else {
			freq[kmer] = 1
		}
	}

	max := 0

	for _, val := range freq {
		if val > max {
			max = val
		}
	}

	most := make([]string, 0, len(freq))

	for kmer, val := range freq {
		if val == max {
			most = append(most, kmer)
		}
	}

	for _, val := range most {
		fmt.Printf("%s ", val)
	}
	fmt.Println()

}
