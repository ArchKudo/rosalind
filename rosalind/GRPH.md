## Problem

**Given**: A collection of [DNA strings](https://rosalind.info/glossary/dna-string/) in [FASTA format](https://rosalind.info/glossary/fasta-format/) having total length at most 10 [kbp](https://rosalind.info/glossary/kbp/).

**Return**: The adjacency list corresponding to O3. You may return edges in any order.

### Sample Dataset

```
>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG
```


### Sample Output

```
Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_2391 Rosalind_2323
```


## Solution

### Naive O($n^2$) Solution

```python
def parse_fasta(file: str):
    with open(file, "r") as f:
        fasta = f.read()
        pairs = [read.splitlines() for read in fasta.split(">")[1:] ]
        pairs = [[x[0], "".join(x[1:])] for x in pairs]
        print(pairs)
        return pairs

def make_graph(pairs: List[List[str]]):
    matches = []

    for i in range(len(pairs)):
        curr = pairs[i]
        back = curr[1][-3:]
        rest = pairs[:i] + pairs[i+1:]
        match = [[curr[0],x[0]] for x in rest if x[1][:3] == back]
        matches += match

    return matches

def printarr(matches: List[List[str]]):
    for i in matches:
        print(f"{i[0]} {i[1]}")

if __name__ == '__main__':
    pairs = parse_fasta("rosalind_grph.txt")
    matches = make_graph(pairs)
    printarr(matches)
```

### Making the loop faster still O($n^2$)

```python
def make_graph(pairs: List[List[str]]) -> List[List[str]]:
	return [
		a[0], b[0]
		for a in pairs
		for b in pairs
		if a is not b and a[1][-3:] == b[1][:3]
	]
```

### Make the loop go faster with dictionary lookups O(n)

```python
def maek_graph(pairs: List[List[str]]) -> List[List[str]]:
    prefixes = {}

    for prefix_id,seq in pairs:
        prefix = seq[:3]
        if prefix not in prefixes:
            prefixes[prefix] = []
        prefixes[prefix].append(prefix_id)

    matches = []
    for suffix_id,seq in pairs:
        suffix = seq[-3:]
        if suffix in prefixes:
            for prefix_id in prefixes[suffix]:
                if prefix_id != suffix_id:
                    matches.append([suffix_id, prefix_id])

    return matches
```

### Benchmark the two

```python

from timeit import timeit

if __name__ == '__main__':
	pairs = parse_fasta("rosalind_grph.txt")
    t = timeit("make_graph(pairs)", globals=locals(), number=1000)
    print(t)
    t = timeit("maek_graph(pairs)", globals=locals(), number=1000)
    print(t)
```