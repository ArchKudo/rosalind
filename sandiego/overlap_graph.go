package main

import (
	"fmt"
	"os"
	"strings"
)

func OverlapGraph(reads []string) map[string][]string {
	L := len(reads)
	sz := len(reads[0])
	lst := make(map[string][]string, L)

	for i := 0; i < L-1; i++ {
		for j := i + 1; j < L; j++ {
			if reads[i][1:] == reads[j][:sz-1] {
				if _, ok := lst[reads[i]]; ok {
					lst[reads[i]] = append(lst[reads[i]], reads[j])
				} else {
					lst[reads[i]] = []string{reads[j]}
				}
			} else if reads[i][:sz-1] == reads[j][1:] {
				if _, ok := lst[reads[j]]; ok {
					lst[reads[j]] = append(lst[reads[j]], reads[i])
				} else {
					lst[reads[j]] = []string{reads[i]}
				}
			}
		}
	}

	return lst
}

func main() {
	file, err := os.ReadFile("overlap_graph.txt")

	var content string

	if err == nil {
		content = string(file)
	} else {
		content = "ATGCG GCATG CATGC AGGCA GGCAT GGCAC"
	}

	reads := strings.Fields(content)

	lst := OverlapGraph(reads)

	for prefix, suffixes := range lst {
		fmt.Printf("%s:", prefix)
		for _, suffix := range suffixes {
			fmt.Printf(" %s", suffix)
		}
		fmt.Println()
	}
}
