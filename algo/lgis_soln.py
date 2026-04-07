def lis(sz: int, arr: list[int]) -> list[int]:

    if sz == 0:
        return arr

    dp = [1] * sz
    prev = [-1] * sz

    for i in range(sz):
        print(f"i={i}")
        for j in range(i):
            print(f"\tj={j}")
            print(f"\t arr[{j}]={arr[j]} arr[{i}]={arr[i]}")
            print(f"\t dp[j]={dp[j]} dp[j]+1={dp[j] + 1} dp[i]={dp[i]}")
            if arr[j] < arr[i] and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1

                print(f"\t\tdp[{i}]={dp[i]}")
                print(f"\t\tdp[{j}]={dp[j]}")
                print(f"\t\tdp[{j}]+1={dp[j] + 1}")

                prev[i] = j
                print(f"\t\tprev[{i}]={prev[i]}")

    bestend = 0

    for i in range(1, sz):
        if dp[i] > dp[bestend]:
            bestend = i

    result = []
    k = bestend

    print(prev)
    while k != -1:
        print(k)
        result.append(arr[k])
        k = prev[k]
        print(f"\t{k}")

    return result[::-1]


if __name__ == "__main__":
    print(lis(5, [5, 1, 4, 2, 3]))
    print(lis(5, [1, 2, 3, 4, 5]))
