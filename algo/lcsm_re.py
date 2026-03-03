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
            print(f"Adding sequence with uid: {idx} and value: {seq}")
            self.add_seq(seq, idx)

        print("Calculating LCSM")
        self.calculate()

        print(f"Found lcsm of length {self.maxl} and value {self.lcsm}")



    def add_seq(self, seq: str, uid: int):
        curr = self.root
        print(f"Current set to root: {curr}")

        self.count += 1

        print(f"Increased count to {self.count}")

        for idx in range(len(seq)):
            print(f"\tWorking on sequence: {seq} from {idx}...{len(seq)}")
            curr = self.root
            print(f"\tReset curr to root: {curr}")
            for char in seq[idx:]:
                print(f"\t\tWorking on {char} from substring: {seq[idx:]}")

                if char not in curr.children:
                    print(f"\t\t\tDidn't find {char} in {curr.children}")
                    curr.children[char] = Node()

                print(f"\t\t{char} already present in {curr.children}")
                curr = curr.children[char]
                curr.covers.add(uid)
                print(f"\t\tCurrent set to {curr}")



    def calculate(self):
        def dfs(node: Node, path: str):
            print(f"Searching through node: {node} with path: {path}")
            print(f"Current {node} covers {node.covers}")
            if len(node.covers) == self.count:
                print(f"\t{node.covers} covers all sequences")
                if len(path) > self.maxl:
                    print(f"\t\tFound current path: {path} longer than {self.maxl}")
                    self.maxl = len(path)
                    self.lcsm = path
                    print(f"\t\tSet LCSM to {path}")

            for char, child in node.children.items():
                print(f"Now running search for {child} and path: {path} + {char}")
                dfs(child, path + char)

        print(f"Count of sequences is {self.count}")
        dfs(self.root, "")


if __name__ == "__main__":

    trie = Trie("GATTACA", "TAGACCA", "ATACA")

    print(trie.lcsm)
