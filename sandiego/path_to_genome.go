package main

import (
	"fmt"
	"os"
	"strings"
)

func GenomePath(fragments []string) string {
	L := len(fragments)
	sz := len(fragments[0])
	genome := fragments[0]

	for i := 1; i < L; i++ {
		genome = genome + string(fragments[i][sz-1])
	}
	return genome
}

func main() {
	file, err := os.ReadFile("path_to_genome.txt")

	var content string

	if err == nil {
		content = string(file)
	} else {
		content = `
			ACCGA CCGAA CGAAG GAAGC AAGCT
		`
	}

	path := strings.Fields(content)

	fmt.Println(GenomePath(path))
}
