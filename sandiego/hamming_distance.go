package main

import (
	"errors"
	"fmt"
	"os"
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

func main() {

	content, err := os.ReadFile(".txt")

	var read, bead string
	if err == nil {
		arr := strings.Split(string(content), "\n")
		read = arr[0]
		bead = arr[1]
	} else {
		read = "CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT"
		bead = "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"
	}

	dist, _ := Hamming(read, bead)
	fmt.Println(dist)
}
