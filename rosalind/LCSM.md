
## Problem

**Given**: A collection of k (k≤100) [DNA strings](https://rosalind.info/glossary/dna-string/) of length at most 1 [kbp](https://rosalind.info/glossary/kbp/) each in [FASTA format](https://rosalind.info/glossary/fasta-format/).

**Return**: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)

### Sample Dataset

```
>Rosalind_1
GATTACA
>Rosalind_2
TAGACCA
>Rosalind_3
ATACA
```

### Sample Output

`AC`

## Initial Thoughts

- Work with smallest string
- Calculate all kmers
	- ABC
	- A B C
	- AB BC
- Search from the largest kmer and early exit on correct match

## Solution

### Bruteforce with extra steps

```python
from typing import List,Set

def parse_fasta(file: str) -> List[str] :
    """
    From previous problem GRPH, which returns both ID / sequence
    The function is changed to return just the sequences,
    making it easier for downstream functions.
    """
    with open(file, "r") as f:
        fasta = f.read()
        pairs = [read.splitlines() for read in fasta.split(">")[1:] ]
        seqs = ["".join(x[1:]) for x in pairs]
        print(seqs)
        return seqs


def kmers(seq: str) -> List[str]:
    """
    Make it easier to visualize the steps for the loops:
    seq[0:1] seq[1:2] seq[2:3] seq[3:4]
    seq[0:2] seq[1:3] seq[2:4]
    seq[0:3] seq[1:4]
    seq[0:4]

	Transpose 
    seq[0:1] seq[0:2] seq[0:3] seq[0:4]
    seq[1:2] seq[1:3] seq[1:4]
    seq[2:3] seq[2:4]
    seq[3:4]

    i -> 0..3 range(len(str))

    k -> 1..4 range(1, len(str) + 1, 1)

    """

    kmers = []

    for i in range(len(seq)):
        for k in range(i + 1, len(seq) + 1, 1):
            kmers.append(seq[i:k])

    return kmers


def kemrs(seq: str) -> Set[str]:
	"""
	Simplify above function
	Remove duplicates
	Sort with longest kmer first
	"""
    return sorted(list({
        seq[i:k]
        for i in range(len(seq))
        for k in range(i + 1, len(seq) + 1, 1)
        }), key=len, reverse=True)


def lcsm(seqs: List[str]) -> str:
    shortest = sorted(seqs, key=len)[0]
    print(shortest)

    kmers = kemrs(shortest)
    print(kmers)

    for kmer in kmers:
        print(kmer)
        res = all([kmer in seq for seq in seqs])
        print(res)

        if res:
            return kmer



if __name__ == "__main__":
    seqs = parse_fasta("rosalind_lcsm.txt")
    # print(kmers("ABCDEFGH"))
    # print(kemrs("ABCD"))

    print(lcsm(seqs))
```


### Simpler brute force

Instead of first calculating all kmers and then finding the solution do both in a single loop.
Also jog my memory of combination of slices by writing down what's expected and printing the for loops until I got it right.
Instead of starting for longest slice and 

```python
def lcsmbf(seqs: List[str]) -> str:
    shortest = min(seqs, key = len)
    
    """
    ABCD
    A [0:1]
    AB [0:2]
    ABC [0:3]
    ABCD [0:4]
    B [1:2]
    BC [1:3]
    BCD [1:4]
    C [2:3]
    CD [2:4]
    D [3:4]
    """
    
    longest = ""

    for i in range(len(shortest)):
        for j in range(i+1, len(shortest)):
            print(f"[{i}:{j}]")
            print(f"{shortest[i:j]}")
            found = all(shortest[i:j] in seq for seq in seqs)
            if found:
                candidate = shortest[i:j]
                print(candidate)
            if len(candidate) > len(longest):
                longest = candidate
                print(longest)


    print(f"Longest is {longest}")

    return longest
```

### A quicker solution using binary search on lengths

Here the binary search is on the length and not the actual substring
So first check at midpoint:
That is all substring of length n/2

If no match is found here it means nothing will be found in longer lengths as well,
reduce the upper bound by to n/2 - 1

If a match is found no point checking lower bound as the longest will be shorter than whatever is the current on, increase the lower bound to n/2 + 1

Continue until upper bound is greater than lower bound

Also do an early break so we only check the first condition which matches for a length as checking the rest for matches will have no effect on longest.

```python
def lcsmbs(seqs: List[str]) -> str:
    shortest = min(seqs, key=len)
    
    sl: int = 1 # shortest possible length of substring
    ll: int = len(shortest) # longest possible length of substring

    best: str = "" # longest common substring


    while sl <= ll:
        
        found = None
        mid = (sl + ll) // 2

        # Prevent going out of bounds
        span = range(len(shortest) - mid + 1)

        for i in span:
            substr = shortest[i: i + mid]
            if all(substr in seq for seq in seqs):
                found = substr
                break

        if found:
            best = found
            sl = mid + 1 # Exclude lengths shorter than current
        else:
            ll = mid - 1 # Exclude lengths longer than current
            
    
    return best
```

### Other solutions

- Also tried using a Suffix Trie instead of a Suffix Tree but that exploded in memory
- Use rolling hashes somehow to make the comparisions faster