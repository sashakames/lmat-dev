import sys

outf = open(sys.argv[2], "w")

next = 0

count = 0


for line in open(sys.argv[1]):
    
    if next == 1:
        outf.write(">11676-line-" + str(count) + "\n")
        outf.write(line)
        next = 0
    if (line.strip() == ">11676"):
        next=1

        print count
    count += 1


                   
