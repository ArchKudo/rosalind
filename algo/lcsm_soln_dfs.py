class Node:
    def __init__(self):
        self.children = {}
        self.indexes = set()

class STree:
    def __init__(self):
        self.root = Node()
        
    def add_suffix(self, suffix, idx):
        node = self.root
        for char in suffix:
            if char not in node.children:
                node.children[char] = Node()
            node = node.children[char]
            node.indexes.add(idx)

    def add_string(self, string, idx):
        for i in range(len(string)):
            self.add_suffix(string[i:], idx)

    def lcsm(self, nstrs):
        self.max_len = 0
        self.result = ""

        def dfs(node, path):
            if len(node.indexes) == nstrs:
                if len(path) > self.max_len:
                    self.max_len = len(path)
                    self.result = path

            for char, child in node.children.items():
                dfs(child, path + char)

        dfs(self.root, "")

        return self.result


def lcsm(strings):
    tree = STree()
    for idx, string in enumerate(strings):
        tree.add_string(string, idx)

    return tree.lcsm(len(strings))

if __name__ == '__main__':
    print(lcsm(["ACAG" , "ACAGGGGA", "GGGGACATTT"]))

