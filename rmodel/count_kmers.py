import sys

f = open(sys.argv[1])


arr = [0] * 21000000

for line in f:
    if line[0] != '#':

        parts = line.split()
        
        count = int(parts[3])

        for i in range(count):
            
            idx = int(parts[4+i*2])
            
            arr[idx] += 1

f.close()
outf = open(sys.argv[1] + ".part.kcnt", "w")        

for i in range(len(arr)):

    if arr[i] > 0:
        outf.write( str(i) + " " + str(arr[i]) + "\n")

outf.close()


