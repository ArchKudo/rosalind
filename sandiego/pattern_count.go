package main

import (
	"fmt"
)

func main() {
	hay, pin := "GACGATATACGACGATA", "ATA"

	count := 0

	for i := 0; i < len(hay)-len(pin)+1; i++ {
		if hay[i:i+len(pin)] == pin {
			count++
		}
	}

	fmt.Println(count)
}
