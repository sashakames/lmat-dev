import sys

f = open(sys.argv[1])

gid_max = 0
tid_max = 0

for line in f:

    
    parts = line.split()

    gid = int(parts[0])
    tid = int(parts[1].strip())


    if (gid > gid_max):
        gid_max = gid

    if(tid > tid_max):
        tid_max = tid

f.close()

print gid_max
print tid_max
