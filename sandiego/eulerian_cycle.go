package main

import (
	"fmt"
	"os"
	"strings"
)

func Cycle(graph map[string][]string, start string) []string {
	i := start
	path := []string{i}

	for {
		stack, ok := graph[i]

		if !ok {
			break
		}

		last := len(stack) - 1
		path = append(path, stack[last])
		graph[i] = stack[:last]

		if len(stack[:last]) == 0 {
			delete(graph, i)
		} else {
			graph[i] = stack[:last]
		}

		i = stack[last]

	}

	return path
}

func CopyMap(graph map[string][]string) map[string][]string {
	copymap := make(map[string][]string, len(graph))

	for k, v := range graph {
		if v == nil {
			continue
		}
		copymap[k] = v
	}

	return copymap
}

func EulerianCycle(graph map[string][]string) []string {
	var start string
	var path []string

	for k, _ := range graph {
		start = k
		break
	}

	for {
		copymap := CopyMap(graph)
		cycle := Cycle(copymap, start)

		fmt.Println(cycle)
		if len(copymap) == 0 {
			path = cycle
			break
		}

		for _, v := range copymap {
			if v != nil {
				start = v[len(v) - 1]
				break
			}
		}

	}

	return path

}

func main() {
	file, err := os.ReadFile("eulerian_cycle.txt")

	var content string

	if err == nil {
		content = string(file)
	} else {
		content = `
			0: 3
			1: 0
			2: 1 6
			3: 2
			4: 2
			5: 4
			6: 5 8
			7: 9
			8: 7
			9: 6
		`
	}

	edges := strings.Split(strings.TrimSpace(content), "\n")
	graph := make(map[string][]string, len(edges))

	for _, edge := range edges {
		lhs, rhs, _ := strings.Cut(strings.TrimSpace(edge), ":")
		graph[lhs] = strings.Fields(rhs)
	}
	
	fmt.Println(graph)
	cycle := EulerianCycle(graph)

	for _, v := range cycle {
		fmt.Printf("%s ", v)
	}
	fmt.Println()

}
