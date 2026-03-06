from functools import reduce
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight


def prtm(pepper: str):

    return molecular_weight(Seq(pepper), seq_type="protein")


def prtm2(pep: str):

    massmap = {
        "A": 71.03711,
        "C": 103.00919,
        "D": 115.02694,
        "E": 129.04259,
        "F": 147.06841,
        "G": 57.02146,
        "H": 137.05891,
        "I": 113.08406,
        "K": 128.09496,
        "L": 113.08406,
        "M": 131.04049,
        "N": 114.04293,
        "P": 97.05276,
        "Q": 128.05858,
        "R": 156.10111,
        "S": 87.03203,
        "T": 101.04768,
        "V": 99.06841,
        "W": 186.07931,
        "Y": 163.06333,
    }

    return round(reduce(lambda x, y: x + massmap[y], pep, initial=0), 3)


if __name__ == "__main__":
    with open("rosalind_prtm.txt", "r") as f:
        pep = f.read().strip()

    print(prtm2(pep))
