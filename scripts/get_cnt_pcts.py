import sys

arr = []

for line in sys.stdin :



    parts = line.strip().split(".")


    if (parts[-1] == "rank_cnt_out"): 

        print (line.strip())

        f = open(line.strip())

        for ll in f:
            
            if (ll[0:3] == "tp=" or ll[0:3] == "fn="):

                print ll.strip()               
            
        f.close()       
    
