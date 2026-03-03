from typing import Dict, List
import re
import urllib.request


def fetch_accession(accession_id: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{accession_id}.fasta"

    return "".join(
                urllib.request.urlopen(url)
                .read().decode("utf-8").splitlines()[1:])


def parse_accessions(filename: str) -> List[str]:
    with open(filename, "r") as f:
        pttrn = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]{1,2}"

        txt = f.read()
        return zip(txt.split("\n"),  re.findall(pttrn, txt))

def mprt(seq: str) -> List[int]:
    motif = re.compile(r"(?=(N[^P][ST][^P]))")

    return [str(match.start(0) + 1) for match in re.finditer(motif, seq)]


def mprt2(seq: str) -> List[int]:
    res = []
    for i in range(len(seq) - 4 + 1):
        group = seq[i:i+4]
        if (group[0] == 'N' and group[1] != 'P' and group[2] in ["S", "T"] and group[3] != 'P'):
            res.append(f"{i + 1}")

    return res



if __name__ == "__main__":


    acc = parse_accessions("rosalind_mprt.txt")

    down = map(lambda x: [x[0], fetch_accession(x[1])], acc)

    res1 = map(lambda x: [x[0], " ".join(mprt(x[1]))], down)
    res2 = map(lambda x: [x[0], " ".join(mprt2(x[1]))], down)

    for a,b in res1:
        if b:
            print(a)
            print(b)

    for a,b in res2:
        if b:
            print(a)
            print(b)
