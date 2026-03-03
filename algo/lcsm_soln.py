from typing import List,Set, Dict

def parse_fasta(file: str) -> List[str] :
    """
    Change to return just the sequences
    """

    with open(file, "r") as f:
        fasta = f.read()
        pairs = [read.splitlines() for read in fasta.split(">")[1:] ]
        seqs = ["".join(x[1:]) for x in pairs]
        # print(seqs)
        return seqs


def kmers(seq: str) -> List[str]:
    """
    seq[0:1] seq[1:2] seq[2:3] seq[3:4]
    seq[0:2] seq[1:3] seq[2:4]
    seq[0:3] seq[1:4]
    seq[0:4]

    seq[0:1] seq[0:2] seq[0:3] seq[0:4]
    seq[1:2] seq[1:3] seq[1:4]
    seq[2:3] seq[2:4]
    seq[3:4]

    i -> 0..3 range(len(str))

    k -> 1..4 range(1, len(str) + 1, 1)

    print([seq[i:k] for k in range(i, len(str) + 1, 1) for i in range(len(str))])
    """

    kmers = []

    for i in range(len(seq)):
        for k in range(i + 1, len(seq) + 1, 1):
            kmers.append(seq[i:k])

    return kmers


def kemrs(seq: str) -> Set[str]:
    return sorted(list({
        seq[i:k]
        for i in range(len(seq))
        for k in range(i + 1, len(seq) + 1, 1)
        }), key=len, reverse=True)


def lcsm(seqs: List[str]) -> str:
    shortest = sorted(seqs, key=len)[0]
    # print(shortest)

    kmers = kemrs(shortest)
    # print(kmers)

    for kmer in kmers:
        # print(kmer)
        res = all([kmer in seq for seq in seqs])
        # print(res)

        if res:
            return kmer

class Node:
    def __init__(self):
        self.children: Dict[str, "Node"] = {}
        self.covers = set()

class Trie:

    def __init__(self, *seqs: str):

        self.root = Node()
        self.lcsm = False
        self.maxl = 0
        self.count = 0

        for idx, seq in enumerate(seqs):
            # print(f"Adding sequence with uid: {idx} and value: {seq}")
            self.add_seq(seq, idx)

        # print("Calculating LCSM")
        self.calculate()

        # print(f"Found lcsm of length {self.maxl} and value {self.lcsm}")



    def add_seq(self, seq: str, uid: int):
        curr = self.root
        # print(f"Current set to root: {curr}")

        self.count += 1

        # print(f"Increased count to {self.count}")

        for idx in range(len(seq)):
            # print(f"\tWorking on sequence: {seq} from {idx}...{len(seq)}")
            curr = self.root
            # print(f"\tReset curr to root: {curr}")
            for char in seq[idx:]:
                # print(f"\t\tWorking on {char} from substring: {seq[idx:]}")

                if char not in curr.children:
                    # print(f"\t\t\tDidn't find {char} in {curr.children}")
                    curr.children[char] = Node()

                # print(f"\t\t{char} already present in {curr.children}")
                curr = curr.children[char]
                curr.covers.add(uid)
                # print(f"\t\tCurrent set to {curr}")



    def calculate(self):
        def dfs(node: Node, path: str):
            # print(f"Searching through node: {node} with path: {path}")
            # print(f"Current {node} covers {node.covers}")
            if len(node.covers) == self.count:
                # print(f"\t{node.covers} covers all sequences")
                if len(path) > self.maxl:
                    # print(f"\t\tFound current path: {path} longer than {self.maxl}")
                    self.maxl = len(path)
                    self.lcsm = path
                    # print(f"\t\tSet LCSM to {path}")

            for char, child in node.children.items():
                # print(f"Now running search for {child} and path: {path} + {char}")
                dfs(child, path + char)

        # print(f"Count of sequences is {self.count}")
        dfs(self.root, "")

def lcsmbf(seqs: List[str]) -> str:
    shortest = min(seqs, key = len)
    """
    ABCD
    A [0:1]
    AB [0:2]
    ABC [0:3]
    ABCD [0:4]
    B [1:2]
    BC [1:3]
    BCD [1:4]
    C [2:3]
    CD [2:4]
    D [3:4]
    """
    longest = ""

    for i in range(len(shortest)):
        for j in range(i+1, len(shortest)):
            print(f"[{i}:{j}]")
            print(f"{shortest[i:j]}")
            found = all(shortest[i:j] in seq for seq in seqs)
            if found:
                candidate = shortest[i:j]
                print(candidate)
            if len(candidate) > len(longest):
                longest = candidate
                print(longest)


    print(f"Longest is {longest}")

    return longest

def lcsmbs(seqs: List[str]) -> str:
    shortest: str = min(seqs, key=len)
    low = 1
    high = len(shortest)
    best = ""

    """
    ABCD
    A [0:1]
    AB [0:2]
    ABC [0:3]
    ABCD [0:4]
    B [1:2]
    BC [1:3]
    BCD [1:4]
    C [2:3]
    CD [2:4]
    D [3:4]
    """

    while low <= high:
        mid = (low + high) // 2
        print(f"Low: {low} ; High: {high} ; mid: {mid}")
        found = None
        span = range(len(shortest) - mid + 1)
        for start in span:
            print(f"Slice is from {start}: {start + mid}")
            substr = shortest[start: start + mid]
            print(f"Substring is {substr}")
            if all(substr in seq for seq in seqs):
                found = substr
                print(f"found: {found}")
                break

        if found:
            best = found
            low = mid + 1
            print(f"New best: {best} ; New low: {low}; high: {high}")
        else:
            high = mid - 1
            print(f"best: {best} ; low: {low}; new high: {high}")

    return best


if __name__ == "__main__":
    seqs = parse_fasta("rosalind_lcsm.txt")
    # print(kmers("ABCDEFGH"))
    # print(kemrs("ABCD"))

    # print(lcsm(seqs))

    # trie = Trie(*seqs)

    # print(trie.lcsm)

    print(lcsmbs(seqs))
