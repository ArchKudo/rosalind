from itertools import product


def lexf(alphabets: str, sz: int):
    lst = sorted(alphabets)
    return map(lambda x: "".join(x), product(lst, repeat=sz))


if __name__ == "__main__":
    with open("rosalind_lexf.txt", "r") as f:
        vals = f.read().splitlines()
        alphabets = "".join(vals[0].split())
        sz = int(vals[1])

        print(alphabets)
        print(sz)

    comb = lexf(alphabets, sz)

    for i in comb:
        print(i)
