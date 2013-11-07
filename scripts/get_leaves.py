import sys

dict = {}
for line in open(sys.argv[1]):

    parts = line.split()
    
    dict[parts[0]] = parts[1:]


def get_leaves(val):

    if dict[val][0] == '0':
        print val
    else:
        print val
        for n in dict[val][1:-1]:
            get_leaves(n)

get_leaves(sys.argv[2])
    
