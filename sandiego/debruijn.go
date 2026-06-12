package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"
)

func Composition(read string, k int) []string {
	L := len(read)
	kmers := make([]string, L-k+1)

	for i := 0; i <= L-k; i++ {
		kmers[i] = read[i : i+k]
	}

	return kmers
}

func DeBruijn(read string, k int) map[string][]string {
	kmers := Composition(read, k)

	L := len(kmers)

	path := make(map[string][]string, L)

	for i := 0; i < L-1; i++ {

		if _, ok := path[kmers[i]]; ok {
			path[kmers[i]] = append(path[kmers[i]], kmers[i+1])
		} else {
			path[kmers[i]] = []string{kmers[i+1]}
		}
	}

	return path

}

func main() {
	file, err := os.ReadFile("debruijn.txt")

	var content string

	if err == nil {
		content = string(file)
	} else {
		content = `
			4
			AAGATTCTCTAAGA
		`
	}

	lst := strings.Fields(content)
	k, _ := strconv.Atoi(lst[0])
	read := lst[1]

	path := DeBruijn(read, k-1)

	for prefix, suffixes := range path {
		fmt.Printf("%s:", prefix)
		for _, suffix := range suffixes {
			fmt.Printf(" %s", suffix)
		}
		fmt.Println()
	}
}
