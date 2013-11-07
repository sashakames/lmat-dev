import sys

if (len(sys.argv) < 2):
    print "Usage: " + sys.argv[0] + " <tcnt-file> > <output-fn>"
    exit (0)

count =2

for line in open(sys.argv[1]):
    
    parts = line.split()

    
    if (parts[0] == "1"):
        print "1 1"
    else:
        print parts[0], count
        count = count + 1

 

