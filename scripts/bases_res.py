import sys

sum = 0


f  = open(sys.argv[1])



for line in f:

    parts = line.split();

    if parts[0] == "read:":
        sum = sum + len(parts[1]) -1 



print sum

