package main

import (
	"fmt"
	"math/rand"
	"os"
	"slices"
	"strconv"
	"strings"
)

type R struct {
	reads []string
	sz    int
	k     int
	L     int
	rng   *rand.Rand
}

func RandSearch(r *R) []string {
	best := Rmers(r)

	for i := 0; i < 1000; i++ {
		seed := Rmers(r)
		for {
			profile := Profile(seed)

			if Score(seed) < Score(best) {
				best = seed
			}

			motifs := Motifs(r, profile)

			if Score(motifs) >= Score(seed) {
				break
			} else {
				seed = motifs
			}
		}

	}
	return best
}

func Rmers(r *R) []string {
	starts := make([]int, r.L, r.L)

	for i := range starts {
		starts[i] = r.rng.Intn(r.sz - r.k)
	}

	rmers := make([]string, r.L, r.L)

	for i, read := range r.reads {
		rmers[i] = read[starts[i] : starts[i]+r.k]
	}

	return rmers

}

func Profile(kmers []string) map[byte][]float64 {
	r := len(kmers)
	c := len(kmers[0])
	l := float64(r + c)
	profile := map[byte][]float64{
		'A': slices.Repeat([]float64{1.0 / l}, c),
		'C': slices.Repeat([]float64{1.0 / l}, c),
		'G': slices.Repeat([]float64{1.0 / l}, c),
		'T': slices.Repeat([]float64{1.0 / l}, c),
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

func Motifs(r *R,
	profile map[byte][]float64) []string {

	memrs := make([]string, r.L)

	for l, read := range r.reads {
		max := -1.0
		for i := 0; i <= r.sz-r.k; i++ {
			kmer := read[i : i+r.k]
			prob := 1.0
			for j, base := range kmer {
				switch base {
				case 'A':
					prob *= profile['A'][j]
				case 'C':
					prob *= profile['C'][j]
				case 'G':
					prob *= profile['G'][j]
				case 'T':
					prob *= profile['T'][j]
				}
			}
			if prob > max {
				max = prob
				memrs[l] = kmer
			}
		}
	}
	return memrs
}

func Score(kmers []string) int {
	r := len(kmers)
	c := len(kmers[0])
	most := make([]string, c)
	score := make([]int, c)
	tot := 0
	for i := 0; i < c; i++ {
		counts := make(map[byte]int, 4)
		for j := 0; j < r; j++ {
			base := kmers[j][i]
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
	file, err := os.ReadFile("ranid_search.txt")

	// content = string(file)

	if err == nil {
		content = string(file)
	} else {

		content = `
			8 5
			CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA
		`

	}
	lst := strings.Fields(content)
	k, _ := strconv.Atoi(lst[0])
	L, _ := strconv.Atoi(lst[1])
	reads := lst[2:]
	sz := len(lst[2])

	r := R{
		reads,
		sz, k, L,
		rand.New(rand.NewSource(42)),
	}

	fmt.Println(RandSearch(&r))

	content = `
		AAGCCAAA
		AATCCTGG
		GCTACTTG
		ATGTTTTG
	`

	reads = strings.Fields(content)

	content = `
		CCA
		CCT
		CTT
		TTG
	`
	kmers := strings.Fields(content)

	fmt.Println(kmers)
	fmt.Println(reads)

	for k, v := range Profile(kmers) {
		fmt.Printf("%c %.2f", k, v)
		fmt.Println()
	}

	r = R{reads, len(reads), 3, 4, rand.New(rand.NewSource(42))}

	fmt.Println(Motifs(&r, Profile(kmers)))
}
