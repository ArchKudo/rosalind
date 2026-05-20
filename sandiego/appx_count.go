package main

import (
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func Hamming(a string, b string) (int, error) {

	arra := strings.Split(a, "")
	brra := strings.Split(b, "")

	if len(arra) != len(brra) {
		return -1, errors.New("Read lengths different!")
	}

	mismatches := 0

	for i, vala := range arra {
		bala := brra[i]

		if vala != bala {
			mismatches += 1
		}
	}

	return mismatches, nil
}

func fuzzy_match(read string, genome string, m int) []int {

	k := len(read)
	L := len(genome)

	pos := make([]int, 0)

	for i := 0; i <= L-k; i++ {

		kmer := genome[i : i+k]
		dist, _ := Hamming(read, kmer)

		if dist <= m {
			pos = append(pos, i)
		}
	}

	return pos
}

func main() {

	content, err := os.ReadFile("appx_count.txt")

	var read, genome string
	var m int

	if err == nil {
		arr := strings.Split(string(content), "\n")
		read = arr[0]
		genome = arr[1]
		m, _ = strconv.Atoi(arr[2])
	} else {
		read = "AAAAA"
		genome = "AACAAGCTGATAAACATTTAAAGAG"
		m = 2
	}
	
	pos := fuzzy_match(read, genome, m)
	fmt.Println(len(pos))

}
