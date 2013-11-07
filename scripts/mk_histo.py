import sys

buckets = [0] * 513

factor = float(sys.argv[1])

offset = 1 - (1/factor)

for line in sys.stdin:

    val = int(line.strip())

    idx = val / factor + offset


for k in buckets:
    print k


    


     
