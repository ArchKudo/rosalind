## Problem

**Given**: A DNA string `s` of length at most 100bp and an array `A` containing at most 20 numbers between 0 and 1.

**Return**: An array B having the same length as `A` in which `B[k]` represents the common logarithm of the probability that a random string constructed with the GC-content found in `A[k]` will match s exactly.

## Concepts

### GC content

If GC content of a read is $x$ for  $x \in [0,1]$ ,
Then probability of a position having the base:

`G/C` is $x/2$
`A/T` is $(1-x)/2$

## Solution

```python
from math import log10
from typing import List
from collections import Counter


def prob(s: str, a: List[int]) -> List[int]:
    counts = Counter(s)
    return [
        (counts["A"] + counts["T"]) * log10((1 - gc_prob) / 2)
        + (counts["G"] + counts["C"]) * log10(gc_prob / 2)
        for gc_prob in a
    ]
```