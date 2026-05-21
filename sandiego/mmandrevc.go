package main

import (
	"fmt"
	mapset "github.com/deckarep/golang-set/v2"
	"os"
	"strconv"
	"strings"
)

func RevC(t string) string {
	x := map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

	rc := make([]byte, 0, len(t))

	for i := range t {
		rc = append(rc, x[t[len(t)-i-1]])
	}

	s := string(rc)

	return s
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

func Neighbours(read string, d uint) mapset.Set[string] {
	bases := mapset.NewSet[string]("A", "C", "G", "T")
	neighbours := mapset.NewSet[string]()
	if d == 0 {
		return neighbours
	}

	if len(read) == 1 {
		return bases
	}

	suffix := Neighbours(read[1:], d)

	si := suffix.Iterator()

	for elem := range si.C {
		if Hamming(elem, read[1:]) < d {
			bi := bases.Iterator()
			for base := range bi.C {
				neighbours.Add(base + elem)
			}
		} else {
			neighbours.Add(string(read[0]) + elem)
		}
	}

	return neighbours
}

func MaxMap(freqs map[string]int) []string {

	mx := 0

	for _, count := range freqs {
		if count > mx {
			mx = count
		}
	}

	most := make([]string, 0)

	for kmer, count := range freqs {
		if count == mx {
			most = append(most, kmer)
		}
	}

	return most
}

func Most(ref string, k uint, d uint) []string {
	freqs := make(map[string]int, 0)

	L := uint(len(ref))

	fmt.Println(L)

	for i := uint(0); i < L-k; i++ {

		nhs := Neighbours(ref[i:i+k], d)

		ni := nhs.Iterator()

		for e := range ni.C {

			if _, ok := freqs[e]; ok {
				freqs[e] += 1
				freqs[RevC(e)] += 1
			} else {
				freqs[e] = 1
				freqs[RevC(e)] += 1
			}
		}
	}

	most := MaxMap(freqs)

	return most
}

func main() {

	// fmt.Println(Hamming("ACG", "ATG"))

	// fmt.Println(Neighbours("ACG", 1))

	// fmt.Println(Most("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1))

	content, err := os.ReadFile("mmrevc.txt")

	var read string
	var k, d uint

	if err == nil {
		arr := strings.Split(string(content), "\n")
		read = arr[0]
		lst := strings.Split(arr[1], " ")

		ki, _ := strconv.Atoi(lst[0])
		k = uint(ki)

		di, _ := strconv.Atoi(lst[1])
		d = uint(di)

	} else {
		read = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
		k = 4
		d = 1
	}

	fmt.Println(Most(read, k, d))
}
