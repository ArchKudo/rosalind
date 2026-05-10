package main

import (
	"fmt"
	"os"
	"strings"
)

func main() {

	content, err := os.ReadFile("reverse_complement.txt")

	var template string

	if err == nil {
		template = string(content)
	} else {
		template = "AAAACCCGGT"
	}

	trans := map[string]string{"A": "T", "T": "A", "C": "G", "G": "C"}

	revcarr := make([]string, len(template))

	for idx := range template {
		comp := trans[string(template[len(template)-idx-1])]

		revcarr = append(revcarr, comp)

	}

	fmt.Println(strings.Join(revcarr, ""))
}
