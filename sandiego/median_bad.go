package main

import (
	"fmt"
	. "github.com/deckarep/golang-set/v2"
	"math"
	"os"
	"strconv"
	"strings"
)

func Kmers(k uint) Set[string] {
	bases := NewSet[string]("A", "C", "G", "T")
	kmers := NewSet[string]()

	if k == 0 {
		return kmers
	}

	if k == 1 {
		return bases
	}

	suffixes := Kmers(k - 1)
	for base := range bases.Iter() {
		for suffix := range suffixes.Iter() {
			kmers.Add(base + suffix)
		}
	}

	return kmers
}

func Hamming(a string, b string) uint {
	d := uint(0)

	for i := range a {
		if a[i] != b[i] {
			d++
		}
	}
	return d
}

func D(kmer string, read string) uint {
	k := len(kmer)
	L := len(read)
	min := uint(math.MaxUint)

	for i := 0; i <= L-k; i++ {
		curr := Hamming(kmer, read[i:i+k])
		if curr < min {
			min = curr
		}
	}

	return min
}

func ΣD(kmer string, reads []string) uint {
	d := uint(0)

	for _, read := range reads {
		d += D(kmer, read)
	}

	return d
}

// May have multiple answers eventhough we just return 1
func Memer(g []string, k uint) string {
	var median string
	min := uint(math.MaxUint)

	for kmer := range Kmers(k).Iter() {
		curr := ΣD(kmer, g)
		if curr < min {
			min = curr
			median = kmer
		}
	}

	return median
}

func main() {
	k := uint(3)
	g := strings.Split("AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTTCGGGACAG", " ")

	content, err := os.ReadFile("median_bad.txt")
	if err == nil {
		lst := strings.Split(string(content), "\n")

		ki, _ := strconv.Atoi(lst[0])
		k = uint(ki)

		g = strings.Split(lst[1], " ")
	}

	fmt.Println(Memer(g, k))

	// fmt.Println(Kmers(0))
	// fmt.Println(Kmers(1))
	// fmt.Println(Kmers(2).Cardinality())
	// fmt.Println(Kmers(3).Cardinality())

	// kmer := "AAA"
	// g := []string{"TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"}
	// fmt.Println(ΣD(kmer, g))

	// k := uint(3)
	// g = strings.Split("AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTTCGGGACAG", " ")

	// fmt.Println(Memer(g, k))
}
