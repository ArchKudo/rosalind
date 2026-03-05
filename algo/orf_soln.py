from typing import Dict, List

from Bio.Data import CodonTable
from Bio.Seq import Seq

table = CodonTable.unambiguous_dna_by_name["Standard"]

# start = table.start_codons
start = ["ATG"]
stop = table.stop_codons


codons = table.forward_table


def translate(dna: str):
    return "".join([codons[dna[i : i + 3]] for i in range(0, len(dna), 3)])


def orf(dna: str, codons: Dict[str, str], start: List[str], stop: List[str]):

    proteins = []
    n = len(dna)

    start_ptr = 0

    while start_ptr <= n - 6:
        if dna[start_ptr : start_ptr + 3] in start:
            for ptr in range(start_ptr + 3, n, 3):
                if dna[ptr : ptr + 3] in stop:
                    proteins.append(translate(dna[start_ptr:ptr]))
                    break

        start_ptr += 1

    return proteins


def orfs(template: str, codons: Dict[str, str], start: List[str], stop: List[str]):
    dna = Seq(template)
    cdna = dna.reverse_complement()

    return set(orf(str(dna), codons, start, stop) + orf(str(cdna), codons, start, stop))


if __name__ == "__main__":
    with open("rosalind_orf.txt", "r") as f:
        fasta = f.read()
        dna = "".join(fasta.splitlines()[1:])

    proteins = orfs(dna, codons, start, stop)

    for protein in proteins:
        print(protein)
