def lis(sz: int, arr: list[int]) -> list[int]:
    """
    The best array holds the longest sequence at that position
    Initially the maximum length at each position is one i.e itself
    """
    best: list[int] = [1] * sz

    """
    Previous holds the index like a linked list
    to the last best sequence
    so we can stitch all to get the longest increasing sequence
    """
    prev: list[int] = [-1] * sz

    for i in range(1, sz):
        for j in range(i):
            bigger = arr[j] < arr[i]

            """
            Assume arr = [1,2,0,3]
            i = 1
            best[1] = best[0] + 1 = 2 _ (I)

            i = 2
            arr[i] < arr[j]
            best[2] = 1 _ (II)
            for both j = 0,1

            i = 3
                j=0 best[3] = best[0] + 1 = 2
                j=1 best[3] = best[1] + 1 = 3 (from I)
                j=2 best[3] = best[2] + 1 = 2 (from II)

            We need to exclude cases where the longest sequence gets shorter

            """
            longer = best[j] + 1 > best[i]

            if bigger and longer:
                best[i] = best[j] + 1
                """
                arr = [5,1,4,2,3]
                prev[best=4] = j = 3
                prev[best=3] = j = 1
                prev[best=1] = j = -1

                stopping after -1 is encountered
                i = 4 -> 3 -> 1 -> stop
                arr[i] = 3 -> 2 -> 1 -> stop
                arr[::-1] = 1,2,3
                """
                prev[i] = j

    farthest = max(range(sz), key=lambda i: best[i])

    seq: list[int] = []

    while farthest != -1:
        seq.append(arr[farthest])
        farthest = prev[farthest]

    return seq[::-1]


if __name__ == "__main__":
    print(lis(5, [5, 1, 4, 2, 3]))
    print(lis(5, [1, 2, 3, 4, 5]))
    print(lis(5, [1, 2, 8, 4, 5]))
    print(lis(4, [1, 2, 0, 3]))
