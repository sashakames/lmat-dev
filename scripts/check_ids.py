# file checks total kmer counts in file from a input list
# argumemt1 counts file, eg tcnt.m9.20.tax_histo
# argument2 list-file contains tids

import sys

dict = {}


for line in open(sys.argv[1]):

    parts = line.split()

    dict[parts[0]] = int(parts[1].strip())


for line in open(sys.argv[2]):

    val = line.strip()

    if val in dict:
        print dict[val]


