package main

import (
	"fmt"
	"math"
	"os"
	"strings"
)

func Skew(read []string) []int {
	skew := make([]int, len(read)+1)

	for i, val := range read {
		switch val {
		case "C":
			skew[i+1] = skew[i] - 1
		case "G":
			skew[i+1] = skew[i] + 1
		default:
			skew[i+1] = skew[i]
		}
	}

	return skew

}

func Minima(skew []int) []int {

	min := math.MaxInt

	for _, val := range skew {
		if val < min {
			min = val
		}
	}

	pos := make([]int, 0)

	for idx, val := range skew {
		if val == min {
			pos = append(pos, idx)
		}
	}

	return pos
}

func main() {
	content, err := os.ReadFile(".txt")

	var read string

	if err == nil {
		read = string(content)
	} else {
		read = string("CATTCCAGTACTTCGATGATGGCGTGAAGA")
	}

	fmt.Println(read)

	bases := strings.Split(read, "")

	fmt.Println(Minima(Skew(bases)))
}
