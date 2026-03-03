## Problem

**Given**: At most 15 UniProt Protein Database access IDs.

**Return**: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of [locations](https://rosalind.info/glossary/location/) in the protein string where the motif can be found.

Sample Dataset

```
A2Z669
B5ZC00
P07204_TRBM_HUMAN
P20840_SAG1_YEAST
```

Sample Output
```
B5ZC00
85 118 142 306 395
P07204_TRBM_HUMAN
47 115 116 382 409
P20840_SAG1_YEAST
79 109 135 248 306 348 364 402 485 501 614
```

## Solution

Wasted entire day trying to use regex but it would give me the wrong answer

The one that works:

```python
def mprt2(seq: str) -> List[int]:
    res = []
    for i in range(len(seq) - 4 + 1):
        group = seq[i:i+4]
        if (group[0] == 'N' and 
	        group[1] != 'P' and 
	        group[2] in "ST" and 
	        group[3] != 'P'):
            res.append(f"{i + 1}")

    return res
```

With helpers:

```python
from typing import List
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
        
if __name__ == "__main__":
    acc = parse_accessions("rosalind_mprt.txt")
    down = map(lambda x: [x[0], fetch_accession(x[1])], acc)
    res2 = map(lambda x: [x[0], " ".join(mprt2(x[1]))], down)
	for a,b in res2:
	if b:
		print(a)
		print(b)
```