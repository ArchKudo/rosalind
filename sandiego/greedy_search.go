package main

import (
	"fmt"
	"strconv"
	"strings"
	"os"
)

/*
 * Inputs
 * reads - Collection of reads
 * k - Length of kmer
 * L - len(reads)
 */
func GreedySearch(reads []string, k int, L int) []string {

	// Assign the first kmer from each read in reads to best
	best := make([]string, L)

	for i := 0; i < L; i++ {
		best[i] = reads[i][0:k]
	}

	// fmt.Println(best)

	// fmt.Println(Score(best))

	// Initialize profile of the kmers motifs
	

	// fmt.Println(profile)

	// fmt.Println(reads[0])
	// fmt.Println(Kmers(reads[0], k))

	for _, motif := range Kmers(reads[0], k) {
		// fmt.Println(motif)
		// Create candidate for best
		motifs := make([]string, 1, L)
		// And add kmer from first reads string into it
		motifs[0] = motif
		fmt.Println(motifs)

		Score(motifs)

		for i := 1; i < L; i++ {
			profile := Profile(motifs)
			fmt.Printf("%.1f\n", profile)
			motifs = append(motifs, MostProbable(reads[i], k, profile))
			fmt.Println(motifs, profile)
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
		kmers[i] = read[i:i+k]
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
			// fmt.Println(string(kmers[j][i]))
			switch kmers[j][i] {
			case 'A':
				// fmt.Println(profile['A'][i])
				profile['A'][i] += 1.0 / float64(r)
				// fmt.Println("v", i, profile['A'][i])
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

			// fmt.Println(prob)
		}

		// fmt.Println(kmer, prob)

		// fmt.Println(kmer)
		// fmt.Println(prob)

		// fmt.Println(prob, max)

		if prob > max {
			max = prob
			memr = kmer
			// fmt.Println(max)
			// fmt.Println(kmer)
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

	fmt.Println(r,c)

	most := make([]string, c) // dbg
	score := make([]int, c)
	tot := 0

	// Iterate through transpose
	for i := 0; i < c; i++ {
		// Count for each base in all rows of a column
		counts := make(map[byte]int, 4)

		// Counts map
		for j := 0; j < r; j++ {
			base := reads[j][i]
			if _, ok := counts[base]; ok {
				counts[base] += 1
			} else {
				counts[base] = 1
			}
		}

		max := 0

		// Find max count
		for base, count := range counts {
			if count > max {
				max = count
				most[i] = string(base) // dbg
				score[i] = r - count
			}
		}

	}


	// fmt.Println(score)

	// Sum of scores
	for _, v := range score {
		tot += v
	}

	return tot
}

func main() {

	var content string

	file, err := os.ReadFile("greedy_searcih.txt")

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

	fmt.Println(reads[0:3], k, L)
	fmt.Println(GreedySearch(reads, k, L))

	// reads := []string{
	// 	"TCGGGGGTTTTT",
	// 	"CCGGTGACTTAC",
	// 	"ACGGGGATTTTC",
	// 	"TTGGGGACTTTT",
	// 	"AAGGGGACTTCC",
	// 	"TTGGGGACTTCC",
	// 	"TCGGGGATTCAT",
	// 	"TCGGGGATTCCT",
	// 	"TAGGGGAACTAC",
	// 	"TCGGGTATAACC",
	// }
	// for _, read := range reads {
	// 	fmt.Println(read)
	// }
	// fmt.Println(Score(reads))

	// read := "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"

	// k := 5

	// profile := map[byte][]float64 {
	// 	'A': {0.2,0.2,0.3,0.2,0.3},
	// 	'C': {0.4,0.3,0.1,0.5,0.1},
	// 	'G': {0.3,0.3,0.5,0.2,0.4},
	// 	'T': {0.1,0.2,0.1,0.1,0.2},
	// }

	// fmt.Println(MostProbable(read, k, profile))


	// k := 4
	// x := Profile([]string{"ACCT"})

	// for _, v := range x {
	// 	fmt.Println(v)
	// }

	
	// n := MostProbable("AGGATCTGTC", k, x)

	
	// x = Profile([]string{"ACCT", n})
	
	// for _, v := range x {
	// 	fmt.Println(v)
	// }

	// n1 := MostProbable("CCGACGTTAG", k, x)

	// x = Profile([]string{"ACCT", n, n1})
	
	// for _, v := range x {
	// 	fmt.Println(v)
	// }
	
	// n2 := MostProbable("CAGCAAGGTG", k, x)

	// fmt.Println(n2)

	// 	x = Profile([]string{"ACCT", n, n1, n2})
		
	// 	for _, v := range x {
	// 		fmt.Println(v)
	// 	}

	
	// 	n3 := MostProbable("CACCTGAGCT", k, x)

	// 	fmt.Println(n3)

	// 		x = Profile([]string{"ACCT", n, n1, n2, n3})
			
	// 		for _, v := range x {
	// 			fmt.Println(v)
	// 		}

	// 	fmt.Println([]string{"ACCT", n, n1, n2, n3})
}
