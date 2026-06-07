package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"
)

func MostProb(read string, k int, profile [][]float64) string {
	var memr string
	max := -1.0
	L := len(read)

	for i := 0; i <= L-k; i++ {
		kmer := read[i : i+k]
		prob := 1.0

		for j, base := range kmer {

			switch base {
			case 'A':
				prob *= profile[0][j]
			case 'C':
				prob *= profile[1][j]
			case 'G':
				prob *= profile[2][j]
			case 'T':
				prob *= profile[3][j]
			}
		}

		if prob > max {
			max = prob
			memr = kmer
		}
	}

	return memr

}

func main() {
	var content string
	file, err := os.ReadFile("kmer_profile.txt")

	if err == nil {
		content = string(file)
	} else {
		content = `
			ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
			5
			0.2 0.2 0.3 0.2 0.3
			0.4 0.3 0.1 0.5 0.1
			0.3 0.3 0.5 0.2 0.4
			0.1 0.2 0.1 0.1 0.2
		`
	}

	lst := strings.Fields(content)
	read := lst[0]
	k, _ := strconv.Atoi(lst[1])

	profile := make([][]float64, 4)

	for i := 0; i < 4; i++ {
		profile[i] = make([]float64, k)
		for j := 0; j < k; j++ {
			ptr := 2 + j + (k * i)
			v, _ := strconv.ParseFloat(lst[ptr], 64)
			profile[i][j] = v
		}
	}

	fmt.Println(MostProb(read, k, profile))
}
