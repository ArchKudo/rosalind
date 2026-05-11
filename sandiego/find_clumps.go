package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"
)

func frequencyTable(window string, k int) map[string]int {
	freq := make(map[string]int, len(window))

	for i := 0; i < len(window)-k+1; i++ {
		kmer := window[i : i+k]

		if _, ok := freq[kmer]; ok {
			freq[kmer] += 1
		} else {
			freq[kmer] = 1
		}
	}

	return freq
}

func main() {
	content, err := os.ReadFile("find_clumps.txt")

	var lst []string
	if err == nil {
		lst = strings.Split(string(content), "\n")
	} else {
		lst = []string{"CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", "5 50 4"}
	}

	genome := lst[0]

	fields := strings.Fields(lst[1])

	arr := make([]int, len(fields))

	for i, v := range fields {
		arr[i], _ = strconv.Atoi(v)
	}
	k := arr[0]
	L := arr[1]
	t := arr[2]

	n := len(genome)

	pattern := make(map[string]struct{})

	for i := 0; i <= n-L; i++ {
		window := genome[i : i+L]
		freqMap := frequencyTable(window, k)
		for kmer, count := range freqMap {
			_, present := pattern[kmer]
			if count >= t && !present {
				pattern[kmer] = struct{}{}
			}
		}
	}

	for k, _ := range pattern {
		fmt.Printf("%s ", k)
	}

}
