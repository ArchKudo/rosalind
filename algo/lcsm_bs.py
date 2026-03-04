from typing import List


def parse_fasta(file: str) -> List[str]:
    """
    Change to return just the sequences
    """

    with open(file, "r") as f:
        fasta = f.read()
        pairs = [read.splitlines() for read in fasta.split(">")[1:]]
        seqs = ["".join(x[1:]) for x in pairs]
        # print(seqs)
        return seqs


def lcsmbs(seqs: List[str]) -> str:
    shortest = min(seqs, key=len)

    sl: int = 1  # shortest possible length of substring
    ll: int = len(shortest)  # longest possible length of substring

    best: str = ""  # longest common substring

    while sl <= ll:
        found = None
        mid = (sl + ll) // 2

        # Prevent going out of bounds
        span = range(len(shortest) - mid + 1)

        for i in span:
            substr = shortest[i : i + mid]
            if all(substr in seq for seq in seqs):
                found = substr
                break

        if found:
            best = found
            sl = mid + 1  # Exclude lengths shorter than current
        else:
            ll = mid - 1  # Exclude lengths longer than current

    return best

    """
    ATACA GATTACA TAGACCA

    low = 1 ; high = 5 ; mid = 3

    [0:3] [1:4] [2:5] - No match found
    
    high = mid - 1 = 2
    low + high // 2 = 1 + 2 // 2 = 3 // 2 = 1

    [start: start + mid]
    
    [0: 1] [1: 2] [2: 3] [3: 4] [4: 5]

    Only 0:1 tested and match found

    Since match found low = mid + 1

    low = 2 ; high = 2; mid = 2 + 2 // 2 = 2

    [start : start + mid] 
    [0:2] [1:3] [2:4] [3:5]

    Only test [0:2] as match found

    low = mid + 1 = 3 > high - stop

    Reduce len to 2

    [0:2] [1:3] [2:4]

    """


if __name__ == "__main__":
    seqs = parse_fasta("lcsm_demo.fa")

    # print(lcsm(seqs))
    # lcsmbf(seqs)
    print(lcsmbs(seqs))
