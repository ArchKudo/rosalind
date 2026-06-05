package main

import (
	"fmt"
	"math"
)

func Entropy(read []byte) float64 {
	counts := make(map[byte]int)
	entropy := float64(0)

	for _, base := range read {
		if _, ok := counts[base]; ok {
			counts[base] += 1
		} else {
			counts[base] = 1
		}
	}

	l := len(read)

	for _, count := range counts {
		p := float64(count) / float64(l)
		e := p * math.Log2(p)
		entropy += e
	}

	return -entropy
}

func Entropies(reads []string) float64 {
	entropy := float64(0)

	for i := range reads[0] {
		col := make([]byte, len(reads))
		for j, read := range reads {
			col[j] = read[i]
		}
		entropy += Entropy(col)
	}

	return entropy
}

func main() {
	reads := []string{
		"TCGGGGGTTTTT",
		"CCGGTGACTTAC",
		"ACGGGGATTTTC",
		"TTGGGGACTTTT",
		"AAGGGGACTTCC",
		"TTGGGGACTTCC",
		"TCGGGGATTCAT",
		"TCGGGGATTCCT",
		"TAGGGGAACTAC",
		"TCGGGTATAACC",
	}

	fmt.Printf("%f\n", Entropies(reads))
}
