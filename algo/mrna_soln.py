from collections import Counter

# from operator import countOf
from math import prod
from typing import Dict

from Bio.Data import CodonTable

table = CodonTable.unambiguous_rna_by_name["Standard"]

aas = table.forward_table.values()

# counts = {aa: countOf(aas, aa) for aa in set(aas)}
counts = Counter(aas)

stops = len(table.stop_codons)


def mrna(protein: str, counts: Dict[str, int], stops: int):
    # lst = list(aa)
    return prod(counts.get(aa, 1) for aa in protein) * stops % 1_000_000


if __name__ == "__main__":
    with open("rosalind_mrna.txt", "r") as f:
        protein = f.read()

    print(mrna(protein, counts, stops))
