from functools import reduce
from typing import List

from Bio.Seq import Seq


def splice(template: str, splices: List[str]):
    return reduce(
        lambda out, intron: out.replace(intron, ""), splices, initial=template
    )


if __name__ == "__main__":
    with open("rosalind_splc.txt", "r") as f:
        fasta = f.read()
        pairs = [read.splitlines() for read in fasta.split(">")[1:]]
        pairs = ["".join(x[1:]) for x in pairs]

    dna = splice(pairs[0], pairs[1:])

    print(dna)

    seq = Seq(dna)

    start = seq.find("ATG")

    print(start)

    print(seq[start:].translate(to_stop=True))
