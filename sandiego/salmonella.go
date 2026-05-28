package main

import (
	"fmt"
	"math"
	"os"
	"strings"

	mapset "github.com/deckarep/golang-set/v2"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/opts"
)

func MinSkew(read string) ([]int, []int) {

	skew := make([]int, len(read)+1)
	pos := make([]int, 0)
	min := math.MaxInt

	for i, val := range strings.Split(read, "") {
		switch val {
		case "C":
			skew[i+1] = skew[i] - 1
		case "G":
			skew[i+1] = skew[i] + 1
		default:
			skew[i+1] = skew[i]
		}
	}

	for _, val := range skew {
		if val < min {
			min = val
		}
	}

	for idx, val := range skew {
		if val == min {
			pos = append(pos, idx)
		}
	}

	return pos, skew[1:]
}

func Plot(read string, skew []int, pos []int, file string) {
	x := strings.Split(read, "")

	y := make([]opts.LineData, 0, len(skew))

	for _, v := range skew {
		y = append(y, opts.LineData{Value: v})
	}

	m := make([]opts.MarkPointNameCoordItem, 0, len(pos))

	for _, idx := range pos {
		m = append(m,
			opts.MarkPointNameCoordItem{
				Name:       "min",
				Coordinate: []interface{}{idx, skew[idx]},
			})
	}

	line := charts.NewLine()

	line.SetXAxis(x).
		AddSeries("skew",
			y,
			charts.WithMarkPointNameCoordItemOpts(m...))

	f, _ := os.Create(file)
	line.Render(f)
}

func RevC(t string) string {
	x := map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	rc := make([]byte, 0, len(t))
	for i := range t {
		rc = append(rc, x[t[len(t)-i-1]])
	}

	return string(rc)
}

func Hamming(a string, b string) uint {
	l := min(len(a), len(b))

	d := uint(0)

	for i := 0; i < l; i++ {
		if a[i] != b[i] {
			d += 1
		}
	}

	return d
}

/**
* Generate all fragments d away from read
* Generate neighbours for read[1:] with hamming < d
* Then Prepend all the bases so hamming = d
* For each of the neighbours of suffix compare hamming with outer
 */
func Neighbours(read string, d uint) mapset.Set[string] {

	bases := mapset.NewSet[string]("A", "C", "G", "T")
	neighbours := mapset.NewSet[string]()

	// Return empty set
	if d == 0 {
		return neighbours
	}

	// e.g. Neighbours(A) = A, C, G, T
	if len(read) == 1 {
		return bases
	}

	suffix := Neighbours(read[1:], d)
	si := suffix.Iterator()

	for e := range si.C {
		if Hamming(e, read[1:]) < d {
			bi := bases.Iterator()
			for base := range bi.C {
				neighbours.Add(base + e)
			}
		} else {
			neighbours.Add(string(read[0]) + e)
		}

		// fmt.Println(e, read, neighbours)
	}

	return neighbours
}

/*
* Given a map of kmers and their occurences
* find the most occuring kmer
 */
func MostFreq(freqs map[string]int) ([]string, int) {
	max := 0

	for _, occ := range freqs {
		if occ > max {
			max = occ
		}
	}

	peaks := make([]string, 0)

	for kmer, occ := range freqs {
		if occ == max {
			peaks = append(peaks, kmer)
		}
	}

	return peaks, max
}

/*
* Find the most occuring (k)-mer
* or identical upto hamming distance < d
* in ref
* for getting potential DnA Box
* for RT to bind near OriC
 */
func DnaABox(ref string, k uint, d uint) ([]string, int) {
	freqs := make(map[string]int, 0)

	L := uint(len(ref))

	for i := uint(0); i < L-k; i++ {

		n := Neighbours(ref[i:i+k], d)

		i := n.Iterator()

		for e := range i.C {
			if _, ok := freqs[e]; ok {
				freqs[e] += 1
				freqs[RevC(e)] += 1
			} else {
				freqs[e] = 1
				freqs[RevC(e)] += 1
			}
		}
	}

	return MostFreq(freqs)
}

func main() {

	content, err := os.ReadFile("salmonella.txt")

	var read string

	if err == nil {

		fasta := strings.Split(string(content), "\n")
		read = strings.Join(fasta[1:], "\n")
	} else {
		read = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
	}

	pos, _ := MinSkew(read)

	fmt.Println(pos)

	// Plot(read, skew, pos, "salmonella_gc_skew.html")

	fmt.Println(RevC(read[:10]))

	fmt.Println(Hamming("GGGCCGTTGGT", "GGACCGTTGAC"))

	fmt.Println(Neighbours("TATA", uint(1)))

	dnaABox, occ := DnaABox(read[pos[0]-500:pos[1]+500],
		uint(9),
		uint(1))

	fmt.Println(dnaABox, occ)
}
