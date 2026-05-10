package main

import (
	"fmt"
	"os"
)

func main() {

	content, _ := os.ReadFile("Vibrio_cholerae.txt")

	pttrn, genome := "CTTGATCAT", string(content)

	var counts []int

	for i := 0; i <= len(genome)-len(pttrn); i++ {
		chunk := genome[i : i+len(pttrn)]

		if chunk == pttrn {
			counts = append(counts, i)
		}
	}

	for _, val := range counts {
		fmt.Printf("%d ", val)
	}

	fmt.Println()
}
