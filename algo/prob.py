from collections import Counter
from math import log10
from typing import List


def prob(s: str, a: List[int]) -> List[int]:
    counts = Counter(s)
    return [
        (counts["A"] + counts["T"]) * log10((1 - gc_prob) / 2)
        + (counts["G"] + counts["C"]) * log10(gc_prob / 2)
        for gc_prob in a
    ]


if __name__ == "__main__":
    s = "CCACATTTCGCTCTATGCCTGAAAACTATCAGCCCTATTGCCTGGCCCGGTTTGTGACTCAGGTTCCAATGCGCGAAGATTGCTG"
    a = list(
        map(
            float,
            "0.063 0.181 0.208 0.263 0.324 0.383 0.443 0.501 0.612 0.636 0.734 0.788 0.843 0.894".split(),
        )
    )
    res = prob(s, a)
    print(*["{0:.3f}".format(r) for r in res])
