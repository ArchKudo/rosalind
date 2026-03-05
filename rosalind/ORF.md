## Problem

**Given**: A [DNA string](https://rosalind.info/glossary/dna-string/) s of length at most 1 [kbp](https://rosalind.info/glossary/kbp/) in [FASTA format](https://rosalind.info/glossary/fasta-format/).

**Return**: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in any order.

### Example

Sample Dataset

```
>Rosalind_99
AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG
```


Sample Output
```
MLLGSFRLIPKETLIQVAGSSPCNLS
M
MGMTPRLGLESLLE
MTPRLGLESLLE
```


### Initial Impressions

1. Get DNA codon table instead of RNA
2. While have not iterated to end of string
	1. if i, i+1, i+2 in start_codon
		1. Iterate i+3, n
			1. if i+3, i+4, i+5 in stop_codon
				1. add to possible proteins

```
XXATGCAATAG

outer loop +1
XXA XAT <ATG>
inner loop + 3:
	CAATAG
```

## Solution

```python
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
```

