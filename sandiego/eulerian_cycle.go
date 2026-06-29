package main

import (
	"fmt"
	"os"
	"strings"
)

func EulerianCycle(graph map[string][]string) []string {
	var start string
	for k := range graph {
		start = k
		break
	}

	stack := []string{start}
	cycle := []string{}

	for len(stack) > 0 {

		// Get the last value (graph-key) from the stack
		v := stack[len(stack) - 1]

		// If there are child nodes to follow
		if len(graph[v]) > 0 {
			// Get the last node
			next := graph[v][len(graph[v]) - 1]
			// Pop the node
			graph[v] = graph[v][:len(graph[v]) - 1]
			// Add the child node to stack
			stack = append(stack, next)
		} else {
			// if graph[v] has no children left add it to cycle
			cycle = append([]string{v}, cycle...)
			// Pop the last item of the stack
			stack = stack[:len(stack) - 1]
		}
	}

	return cycle
}

func main() {
	file, err := os.ReadFile("eulerian_cycl.txt")

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

	cycle := EulerianCycle(graph)

	fmt.Println(cycle)
}
