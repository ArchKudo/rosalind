package main

import (
	"fmt"
	"os"
	"strings"
)

func DeBruijn(kmers []string) map[string][]string {
	L := len(kmers)
	k := len(kmers[0]) - 1

	path := make(map[string][]string, L-k+1)

	for i := 0; i < L; i++ {
		prefix := kmers[i][:k]
		suffix := kmers[i][1:]
		if _, ok := path[prefix]; ok {
			path[prefix] = append(path[prefix], suffix)
		} else {
			path[prefix] = []string{suffix}
		}
	}

	return path
}

func PrintPath(path map[string][]string) {
	for prefix, suffixes := range path {
		fmt.Printf("%s:", prefix)
		for _, suffix := range suffixes {
			fmt.Printf(" %s", suffix)
		}
		fmt.Println()
	}
}

func main() {

	file, err := os.ReadFile("debruijn_kmers.txt")

	var content string

	if err == nil {
		content = string(file)
	} else {
		content = "GAGG CAGG GGGG GGGA CAGG AGGG GGAG"
	}

	kmers := strings.Fields(content)

	path := DeBruijn(kmers)

	PrintPath(path)
}
