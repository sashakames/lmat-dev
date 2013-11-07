from sys import argv


outarr = []

N = int(argv[1])
arr = [0, N-1, N/2, (N/2) -1]


if (len(argv) > 3):
    S1 = int(argv[2])
    S2 = int(argv[3])
else:
    S1 = 0
    S2 = N

for n in range(N/4):
    for i in range(4):
        outarr.append(arr[i])

        if i % 2 == 0:
            arr[i] = arr[i] + 1
        else:
            arr[i] = arr[i] - 1



for i in range(S1, S2):
    print outarr[i]
