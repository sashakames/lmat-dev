import sys

if len(sys.argv) < 2:
    print "acsii_search.py:  utility to brute-force search for kmers in tax_histo ascii output files"
    print "Usage: <list-of-files>  | python ascii_search.py <kmer-list-file>"
    print "<kmer-list-file> contains 64-bit decimal integers"
    print "example:  ls /home/ames4/vbonly/*.db | python ..." 

arr = []

f2 = open(sys.argv[1])

for line in f2:
    arr.append(line.rstrip())
    
f2.close()
count = 0

for ll in sys.stdin:

    f = open(ll.strip())

    for line in f:
        parts = line.split()
        
        for x in arr:

            if (x == parts[0]):
                print line
                count = count + 1
                if (count == len(arr)):
                    exit(0)

    f.close()


