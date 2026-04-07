from itertools import combinations
from typing import Callable


def inc(vec: list[int]) -> bool:
    return all(nxt >= cur for cur, nxt in zip(vec, vec[1:]))


def dec(vec: list[int]) -> bool:
    return all(cur >= nxt for cur, nxt in zip(vec, vec[1:]))


def nxt(mat: list[list[int]], fun: Callable[[list[int]], bool]) -> list[int] | None:
    return next(filter(fun, mat), None)


def lxs(sz: int, vec: list[int], fun: Callable[[list[int]], bool]) -> list[int]:
    if fun(vec):
        return vec
    for i in range(sz - 1, 0, -1):
        if x := nxt(combinations(vec, i), fun):
            return x


def lgis(sz: int, vec: list[int]) -> tuple(list[int], list[int]):
    return (lxs(sz, vec, inc), lxs(sz, vec, dec))


if __name__ == "__main__":
    with open("rosalind_lgis.txt", "r") as f:
        tup = f.read().splitlines()
        sz = int(tup[0])
        lst = list(map(int, tup[1].split(" ")))

    for seq in lgis(sz, lst):
        print(" ".join(str(x) for x in seq))
