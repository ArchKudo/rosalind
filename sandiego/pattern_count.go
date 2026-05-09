package main

import (
	"fmt"
	"log"
	"os"
	"strings"
)

func main() {
	content, err := os.ReadFile("pattern_count.txt")
	if err != nil {
		log.Fatal("err")
	}

	lst := strings.Split(string(content), "\n")
	hay, pin := lst[0], lst[1]

	count := 0

	for i := 0; i < len(hay)-len(pin)+1; i++ {
		if hay[i:i+len(pin)] == pin {
			count++
		}
	}

	fmt.Println(count)
}
