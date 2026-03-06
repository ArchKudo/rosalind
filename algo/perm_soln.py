from itertools import permutations


def fact(n):
    if n == 0 or n == 1:
        return 1

    return n * fact(n - 1)


def perm(n: int):
    print(fact(n))

    for tup in permutations(range(1, n + 1), n):
        for elem in tup:
            print(elem, end=" ")
        print()


if __name__ == "__main__":
    perm(6)
