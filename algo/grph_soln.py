from timeit import timeit
from typing import List


def parse_fasta(file: str) -> List[List[str]]:
    """
    Inputs:
    file: str - Name of fasta file containing sequences

    Output:
    pairs: List[List[str]] - A list containing a list of id,sequence pairs

    e.g.

    >Seq1
    AAA
    BBB

    [["Seq1", "AAABBB"]]

    """
    with open(file, "r") as f:
        fasta = f.read()
        pairs = [read.splitlines() for read in fasta.split(">")[1:]]
        pairs = [[x[0], "".join(x[1:])] for x in pairs]
        return pairs


def make_graph(pairs: List[List[str]]) -> List[List[str]]:
    return [
        [a[0], b[0]]
        for a in pairs
        for b in pairs
        if a is not b and a[1][-3:] == b[1][:3]
    ]


def maek_graph(pairs: List[List[str]]) -> List[List[str]]:
    prefixes = {}

    for prefix_id, seq in pairs:
        prefix = seq[:3]
        if prefix not in prefixes:
            prefixes[prefix] = []
        prefixes[prefix].append(prefix_id)

    matches = []
    for suffix_id, seq in pairs:
        suffix = seq[-3:]
        if suffix in prefixes:
            for prefix_id in prefixes[suffix]:
                if prefix_id != suffix_id:
                    matches.append([suffix_id, prefix_id])

    return matches


def printarr(matches: List[List[str]]):
    for i in matches:
        print(f"{i[0]} {i[1]}")


if __name__ == "__main__":
    pairs = parse_fasta("rosalind_grph.txt")
    # printarr(make_graph(pairs))
    printarr(maek_graph(pairs))

    t = timeit("make_graph(pairs)", globals=locals(), number=1000)
    print(t)
    t = timeit("maek_graph(pairs)", globals=locals(), number=1000)
    print(t)
