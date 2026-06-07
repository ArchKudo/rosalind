package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"
)

/*
 * Inputs
 * reads - Collection of reads
 * k - Length of kmer
 * L - len(reads)
 */
func GreedySearch(reads []string, k int, L int) []string {
	best := make([]string, L)
	for i := 0; i < L; i++ {
		best[i] = reads[i][0:k]
	}
	for _, motif := range Kmers(reads[0], k) {
		motifs := make([]string, 1, L)
		motifs[0] = motif
		Score(motifs)
		for i := 1; i < L; i++ {
			profile := Profile(motifs)
			motifs = append(motifs, MostProbable(reads[i], k, profile))
		}
		if Score(motifs) < Score(best) {
			best = motifs
		}
	}
	return best
}
func Kmers(read string, k int) []string {
	L := len(read)
	kmers := make([]string, L-k+1)
	for i := 0; i <= L-k; i++ {
		kmers[i] = read[i : i+k]
	}
	return kmers
}

/*
 * Update profile given original profile, and a new kmer
 * A profile is a sparse matrix with
 * rows as the 4 bases and
 * columns as the probability of the base
 * appearing at that position across all kmers
 */
func Profile(kmers []string) map[byte][]float64 {
	r := len(kmers)
	c := len(kmers[0])
	profile := map[byte][]float64{
		'A': make([]float64, c),
		'C': make([]float64, c),
		'G': make([]float64, c),
		'T': make([]float64, c),
	}
	for i := 0; i < c; i++ {
		for j := 0; j < r; j++ {
			switch kmers[j][i] {
			case 'A':
				profile['A'][i] += 1.0 / float64(r)
			case 'C':
				profile['C'][i] += 1.0 / float64(r)
			case 'G':
				profile['G'][i] += 1.0 / float64(r)
			case 'T':
				profile['T'][i] += 1.0 / float64(r)
			}
		}
	}
	return profile
}

/*
 * Find kmer of a read which has the highest probability
 * given a profile
 */
func MostProbable(read string,
	k int,
	profile map[byte][]float64) string {
	var memr string
	max := -1.0
	for _, kmer := range Kmers(read, k) {
		prob := 1.0
		for i, base := range kmer {
			switch base {
			case 'A':
				prob *= profile['A'][i]
			case 'C':
				prob *= profile['C'][i]
			case 'G':
				prob *= profile['G'][i]
			case 'T':
				prob *= profile['T'][i]
			}
		}
		if prob > max {
			max = prob
			memr = kmer
		}
	}
	return memr
}

/*
 * Sum of count of occurences of alternate bases for each position given a collection of kmers
 */
func Score(reads []string) int {
	r := len(reads)
	c := len(reads[0])
	most := make([]string, c)
	score := make([]int, c)
	tot := 0
	for i := 0; i < c; i++ {
		counts := make(map[byte]int, 4)
		for j := 0; j < r; j++ {
			base := reads[j][i]
			if _, ok := counts[base]; ok {
				counts[base] += 1
			} else {
				counts[base] = 1
			}
		}
		max := 0
		for base, count := range counts {
			if count > max {
				max = count
				most[i] = string(base)
				score[i] = r - count
			}
		}
	}
	for _, v := range score {
		tot += v
	}
	return tot
}
func main() {
	var content string
	file, err := os.ReadFile("greedy_search.txt")
	if err == nil {
		content = string(file)
	} else {
		content = `
			3 5
			GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG
		`
	}
	lst := strings.Fields(content)
	k, _ := strconv.Atoi(lst[0])
	L, _ := strconv.Atoi(lst[1])
	reads := lst[2:]

	fmt.Println(GreedySearch(reads, k, L))
}
