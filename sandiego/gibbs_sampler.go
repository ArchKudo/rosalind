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
	N     int
	rng   *rand.Rand
}

func GibbsSearch(r *R) []string {
	best := Rmers(r)
	for i := 0; i < 20; i++ {
		seed := Rmers(r)
		for n := 0; n < r.N; n++ {
			i := r.rng.Intn(r.L)
			rest := slices.Clone(seed[:i])
			rest = append(rest, seed[i+1:]...)
			profile := Profile(rest)
			seed[i] = Probs(i, r, profile)
			if Score(seed) < Score(best) {
				fmt.Println(Score(seed), seed)
				best = slices.Clone(seed)
			}
		}
	}
	return best
}

func Probs(idx int, r *R, profile map[byte][]float64) string {
	l := r.sz - r.k + 1
	read := r.reads[idx]
	kmers := make([]string, l)
	for i := 0; i < l; i++ {
		kmers[i] = read[i : i+r.k]
	}
	probs := make([]float64, l)
	for i, kmer := range kmers {
		cum := 1.0
		for j, base := range kmer {
			switch base {
			case 'A':
				cum *= profile['A'][j]
			case 'C':
				cum *= profile['C'][j]
			case 'G':
				cum *= profile['G'][j]
			case 'T':
				cum *= profile['T'][j]
			}
		}
		probs[i] = cum
	}
	for i := 1; i < l; i++ {
		probs[i] += probs[i-1]
	}
	c := r.rng.Float64() * probs[l-1]
	var sel string
	for i, prob := range probs {
		if c <= prob {
			sel = kmers[i]
			return sel
		}
	}
	return sel
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
	file, err := os.ReadFile("gibbs_sampler.txt")
	if err == nil {
		content = string(file)
	} else {
		content = `
			8 5 100
			CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA
			`
	}
	lst := strings.Fields(content)
	k, _ := strconv.Atoi(lst[0])
	L, _ := strconv.Atoi(lst[1])
	N, _ := strconv.Atoi(lst[2])
	reads := lst[3:]
	sz := len(lst[3])
	r := R{
		reads,
		sz, k, L, N,
		rand.New(rand.NewSource(1)),
	}
	fmt.Println(GibbsSearch(&r))
}
