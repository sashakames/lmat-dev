import sys

fname = sys.argv[1]

lenf = int(sys.argv[2])

nreads = int(sys.argv[3])

count = 0

num = lenf / nreads

if (num % 2  == 1):
    num += 1


for line in open(sys.argv[1]):

    if (count < 2):
        print line.strip()
    
    count += 1
    if count == num:
        count = 0
