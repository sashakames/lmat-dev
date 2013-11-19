import sys

global last_tup
global parts

odd = False

last_tup = []

rtotal = 4000000

parts = []



def ret(x):
    global last_tup
    global parts

    if (x > (len(parts) -1)):
            print stderr, "error: parts ", x, ' ', len(parts)
            exit(1)

    if (x > (len(last_tup) -1)):
            print "error: tup ", x , ' ', len(last_tup)
            exit(1)
        

    return float(parts[x]) + float(last_tup[x])

for line in sys.stdin:
    

    parts = line.strip().split('\t')

    
    if (odd):
        odd = False
        
        total = ret(5)

        prec = 1 -(( ret(9) + ret(14) + ret(12)     ) / total)

        recall = ( ret(6) + ret(7) + ret(10) + ret(13)) / rtotal
        nn = (ret(8) + ret(11)) / rtotal
        
        print parts[0],  parts[1],  parts [2],  prec,  recall,  nn

    else:
        odd = True
        last_tup = parts

    
