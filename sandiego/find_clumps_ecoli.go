package main

import (
	"fmt"
	"os"
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
	content, _ := os.ReadFile("E_coli.txt")

	genome := string(content)

	k := 9
	L := 500
	t := 3

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

	fmt.Println(len(pattern))

}
