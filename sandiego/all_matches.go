package main

import (
	"fmt"
	"os"
	"strings"
)

func reverse_complement(template string) string {

	trans := map[byte]byte{
		'A': 'T',
		'T': 'A',
		'C': 'G',
		'G': 'C',
	}

	revc := make([]byte, len(template))

	for idx := range template {
		revc[idx] = trans[template[len(template)-idx-1]]
	}

	return string(revc)

}

func main() {

	content, err := os.ReadFile("all_matches.txt")

	var lst []string

	if err == nil {
		lst = strings.Split(string(content), "\n")
	} else {
		lst = []string{"ATAT", "GATATATGCATATACTT"}
	}

	pttrn, genome := lst[0], lst[1]

	revc := reverse_complement(pttrn)

	var counts []int

	for i := 0; i <= len(genome)-len(pttrn); i++ {
		chunk := genome[i : i+len(pttrn)]

		if chunk == pttrn || chunk == revc {
			counts = append(counts, i)
		}
	}

	for _, val := range counts {
		fmt.Printf("%d ", val)
	}

	fmt.Println()
}
