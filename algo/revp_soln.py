from Bio.Seq import Seq


def revp(dna: str):

    seq53 = Seq(dna)
    rev35 = seq53.reverse_complement()[::-1]

    pals = [
        (i + 1, j)
        for j in range(4, 13)
        for i in range(len(seq53) - j + 1)
        if seq53[i : i + j] == rev35[i : i + j][::-1]
    ]

    return pals


if __name__ == "__main__":
    with open("rosalind_revp.txt", "r") as f:
        seq = "".join(f.read().splitlines()[1:])

    for i, j in revp(seq):
        print(f"{i} {j}")
