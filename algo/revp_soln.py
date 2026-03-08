from typing import List
from Bio.Seq import Seq


def revp(dna: str) -> List[(int, int)]:

    seq = Seq(dna)

    pals = [
        (i + 1, j)
        for j in range(4, 13)
        for i in range(len(seq) - j + 1)
        if seq[i : i + j] == seq[i : i + j].reverse_complement()
    ]

    return pals


if __name__ == "__main__":
    with open("rosalind_revp.txt", "r") as f:
        seq = "".join(f.read().splitlines()[1:])

    for i, j in revp(seq):
        print(f"{i} {j}")
