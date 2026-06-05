package main

import (
	"fmt"
	. "github.com/deckarep/golang-set/v2"
	"os"
	"strconv"
	"strings"
)

var k, d uint
var reads []string

func Hamming(a string, b string) uint {
	d := uint(0)
	l := min(len(a), len(b))

	for i := 0; i < l; i++ {
		if a[i] != b[i] {
			d += 1
		}
	}

	return d
}

func Neighbours(read string) Set[string] {
	bases := NewSet[string]("A", "C", "G", "T")
	neighbours := NewSet[string]()

	if d == 0 {
		return neighbours
	}
	if len(read) <= 1 {
		return bases
	}

	suffix := Neighbours(read[1:])

	for val := range suffix.Iter() {
		if Hamming(val, read[1:]) < d {
			for base := range bases.Iter() {
				neighbours.Add(base + val)
			}
		} else {
			neighbours.Add(string(read[0]) + val)
		}
	}

	return neighbours
}

func Kmers(read string) Set[string] {
	kmers := NewSet[string]()
	l := uint(len(read))

	for i := uint(0); i <= l-k; i++ {
		kmers.Add(read[i : i+k])
	}

	return kmers
}

/*
 * Check if pttrn fuzzy matches all reads in list
 */
func All(pttrn string, reads []string) bool {
	for _, read := range reads {
		// Any of the fuzzy fails return false
		if !Fuzzy(pttrn, read) {
			return false
		}
	}
	return true
}

/*
 * Check if pttrn is at Hamming dance d from any of the k-mer's of read
 * ATA, TATA -> ATA, [TAT, ATA] -> true
 */
func Fuzzy(pttrn string, read string) bool {
	kmers := Kmers(read)
	for kmer := range kmers.Iter() {
		if Hamming(kmer, pttrn) <= d {
			return true
		}
	}
	return false
}

func Motif(reads []string) Set[string] {
	common := NewSet[string]()

	for kmer := range Kmers(reads[0]).Iter() {
		for neighbour := range Neighbours(kmer).Iter() {
			if All(neighbour, reads) {
				common.Add(neighbour)
			}
		}
	}

	return common
}

func main() {
	content, err := os.ReadFile("motif_bad.txt")

	if err == nil {
		lst := strings.Split(string(content), "\n")
		params := strings.Split(lst[0], " ")
		ki, _ := strconv.Atoi(params[0])
		k = uint(ki)
		di, _ := strconv.Atoi(params[1])
		d = uint(di)
		reads = strings.Split(lst[1], " ")
	} else {
		reads = []string{"ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"}
		k = 3
		d = 1
	}

	for motif := range Motif(reads).Iter() {
		fmt.Printf("%s ", motif)
	}
}
