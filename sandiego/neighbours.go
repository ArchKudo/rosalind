package main

import (
	"fmt"
	mapset "github.com/deckarep/golang-set/v2"
	"os"
	"strconv"
	"strings"
)

func Hamming(a string, b string) int {

	d := 0

	l := min(len(a), len(b))

	for i := 0; i < l; i++ {
		if a[i] != b[i] {
			d += 1
		}
	}

	return d
}

func Neighbours(pttrn string, d int) mapset.Set[string] {

	bases := mapset.NewSet[string]("A", "C", "G", "T")

	if d == 0 {
		return mapset.NewSet[string]()
	}

	if len(pttrn) == 1 {
		return bases
	}

	neighbours := mapset.NewSet[string]()

	suffix := Neighbours(pttrn[1:], d)

	for val := range mapset.Elements(suffix) {
		if Hamming(val, pttrn[1:]) < d {
			for base := range mapset.Elements(bases) {
				neighbours.Add(base + val)
			}
		} else {
			neighbours.Add(string(pttrn[0]) + val)
		}
	}

	return neighbours
}

func main() {

	// read := "ACG"
	// d := 1

	content, err := os.ReadFile(".txt")

	var read string
	var d int

	if err == nil {
		x := strings.Split(string(content), "\n")
		read = x[0]
		d, _ = strconv.Atoi(x[1])
	} else {
		read = "ACGT"
		d = 4
	}

	neighbours := Neighbours(read, d)

	for val := range mapset.Elements(neighbours) {
		fmt.Printf("%s ", val)
	}

	fmt.Println(neighbours.Cardinality())

}
