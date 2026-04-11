def overlap(left: str, right: str) -> int:
    minl = len(left) // 2 + 1
    maxl = min(len(left), len(right))

    for s in range(maxl, minl - 1, -1):
        """
        left  = AAAA
        right = AAAA
        
        minl = 4 // 2 + 1 = 3
        maxl = min(4,4) = 4
        
        range(4, 3, -1)

        left[-4:] = left = "AAAA"
        right[:4] = right = "AAAA"

        s = 4

        left =  ATTAGACCTG
        right = ...AGACCTGCCG

        minl = 10 // 2 + 1 = 6
        maxl = min(10,10) = 10

        range(10, 6, -1)

        left[-10:] = ATTAGACCTG
        right[:10] = AGACCTGCCG

        left[-9:] = TTAGACCTG
        right[:9] = AGACCTGCC

        left[-8:] = TAGACCTG
        right[:8] = AGACCTGC

        left[-7:] = AGACCTG
        right[:7] = AGACCTG
        
        s = 7

        1 ATTAGACCTG
        2 ...AGACCTGCCG
        3 ......CCTGCCGGAA
        4 .........GCCGGAATAC
        
        ov(1,2) = 7 ov(2, 1) = 0 ov(3,1) = 0 ov(4,1) = 0
        ov(1,3) = 0 ov(2, 3) = 7 ov(3,2) = 0 ov(4,2) = 0
        ov(1,4) = 0 ov(2, 4) = 0 ov(3,4) = 7 ov(4,3) = 0

        """
        if left[-s:] == right[:s]:
            return s

    return 0


def long(seqs: list[str]) -> str:
    """
    1 ATTAGACCTG
    2 ...AGACCTGCCG
    3 ......CCTGCCGGAA
    4 .........GCCGGAATAC

    ov(1,2) = 7 ov(2, 1) = 0 ov(3,1) = 0 ov(4,1) = 0
    ov(1,3) = 0 ov(2, 3) = 7 ov(3,2) = 0 ov(4,2) = 0
    ov(1,4) = 0 ov(2, 4) = 0 ov(3,4) = 7 ov(4,3) = 0

    rights  = {1: (2, 7), 2: (3, 7), 3: (4, 7)}

    hasleft(1,2) = [0, 1, 0, 0]
    hasleft(2,3) = [0, 1, 1, 0]
    hasleft(3,4) = [0, 1, 1, 1]
    """

    # Dictionary storing the next, overlap information for given index
    rights: dict[int, tuple[int, int]] = {}

    # Find the left-most sequence
    hasleft: list[int] = [False] * len(seqs)

    # O(n * n -1) = O(n^2)
    for lidx, left in enumerate(seqs):
        for ridx, right in enumerate(seqs):
            if lidx == ridx:
                continue

            # O(L^2)
            ov = overlap(left, right)

            if ov:
                hasleft[ridx] = True
                rights[lidx] = (ridx, ov)

            # Tot: O(n^2*L^2)

    # O(n)
    start = next(idx for idx, count in enumerate(hasleft) if not count)
    contig = seqs[start]
    curr = start

    # Foot gun: Here while curr cannot be used and will exit if curr is 0
    # Hence check membership instead
    # O(nL)
    while curr in rights:
        nxt, ov = rights[curr]
        contig += seqs[nxt][ov:]
        curr = nxt

    # O(n^2*L^2) + O(nL) + O(n)
    # where n = len(seqs) where L = len(max, key=len))
    return contig


def parse_fasta(file: str) -> list[str]:
    """
    Inputs:
    file: str - Name of fasta file containing sequences

    Output:
    fragments: List[str] - A list containing sequences

    e.g.

    >Seq1
    AAA
    BBB
    >Seq2
    BBB
    CCC

    ["AAABBB", "BBBCCC"]

    """
    with open(file, "r") as f:
        fasta = f.read()
        pairs = [read.splitlines() for read in fasta.split(">")[1:]]
        return ["".join(x[1:]) for x in pairs]


if __name__ == "__main__":
    seqs = parse_fasta("rosalind_long.txt")

    print(long(seqs))
