import sys

arr = []

for line in sys.stdin :



    parts = line.strip().split(".")


    if (parts[-1] == "check" and parts[-3] == "summary"):

        print (line.strip())

        f = open(line.strip())

        for ll in f:
            
            if (ll[0:6] == "final:"):

                print "param: " + parts[-2] + " "  + ll.strip()               
            
        f.close()       
    
