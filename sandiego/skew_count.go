package main

import (
	"fmt"
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

func main() {

	read := strings.Split("GAGCCACCGCGATA", "")

	fmt.Println(read)
	fmt.Println(Skew(read))

	read = strings.Split("CATGGGCATCGGCCATACGCC", "")
	fmt.Println(read)
	fmt.Println(Skew(read))
}
