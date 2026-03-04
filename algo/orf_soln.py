import re
from typing import Dict, List

from Bio.Data import CodonTable
from Bio.Seq import Seq

table = CodonTable.unambiguous_dna_by_name["Standard"]

start = table.start_codons
stop = table.stop_codons


codons = table.forward_table


def translate(dna: str):
    return "".join([codons[dna[i : i + 3]] for i in range(0, len(dna) - 3, 3)])


def orf(dna: str, codons: Dict[str, str], start: List[str], stop: List[str]):

    seq = Seq(dna)
    rev = seq.reverse_complement()

    pttrn = re.compile(
        rf"""
        (?:{"|".join(start)})
        (?:(?!{"|".join(stop)})[ACGT]{{3}})*
        (?:{"|".join(stop)})
        """,
        re.VERBOSE,
    )

    def trans(seq):
        return [translate(orf) for orf in re.findall(pttrn, str(seq))]

    return set(trans(seq) + trans(rev))


if __name__ == "__main__":
    dna = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"

    print(orf(dna, codons, start, stop))
