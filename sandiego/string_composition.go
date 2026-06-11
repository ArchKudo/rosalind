package main

import (
	"fmt"
	"os"
	"strings"
	"strconv"
)

func Composition(read string, k int) []string {
	L := len(read)
	kmers := make([]string, L - k + 1)

	for i := 0; i <= L - k; i++ {
		kmers[i] = read[i:i+k]
	}

	return kmers
}

func main() {
	file, err := os.ReadFile("string_composition.txt")

	var content string

	if err == nil {
		content = string(file)
	} else {
		content = `
			5
			CAATCCAAC
		`
	}

	lst := strings.Fields(content)
	k, _ := strconv.Atoi(lst[0])
	read := lst[1]

	fmt.Println(Composition(read, k))
}
