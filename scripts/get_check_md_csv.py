import sys


dict = {}


out_fn = ""

vals = range(1,301)


for fn in sys.argv[1:5]:

    

    for i in vals:

        f = open(fn + str(i) + '.check')


        
        for ll in f:
            if (ll[0:6] == "final:"):
                        
                llp = ll.split()

                if (i in dict.keys()):
                    row = dict[i]
                else:
                    row = []
                

                row.append(llp[2])
                row.append(llp[4])
                row.append(llp[6])
                row.append(llp[8])

                dict[i] = row


for fn in sys.stdin:  # metaphlan

    i = -1

    parts = fn.split('.')


    if (len(parts) ==8):

        msv = 0

        if  (len(parts[-3]) > 0):
            msv = 100 * int(parts[-3])

        i = int(msv) + int(parts[-2])
        
        f = open(fn.strip())

        print i

        for ll in f:
            if (ll[0:6] == "final:"):
                
                llp = ll.split()                

                row = dict[i]


                row.append(llp[2])
                row.append(llp[4])
                row.append(llp[6])
                row.append(llp[8])

                dict[i] = row

                break
        f.close()



    
outf = open("all_check.csv" , "w")
    
vv = dict.keys()

vv.sort(key=lambda f: int(f))


for val in vv:



    outf.write(str(val))
        

    row = dict[val]

    for cell in row:
        
        outf.write("," + cell)
            
    outf.write("\n")

outf.close()


        
    





    
