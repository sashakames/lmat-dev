import sys

arr = []

dbs = ["marker_lib_17", "marker_lib", "marker_lib_20"]


dict = {}


out_fn = ""


for line in sys.stdin :

    parts = line.strip().split(".")

    if (out_fn == ""):
        for pp in parts[0:int(sys.argv[1])]:
            out_fn = out_fn + pp + "."
        print out_fn



    db_name = ""

    for pp in parts:

        for dd in dbs:

            if pp == dd:
                    
                db_name = pp
            
        if not db_name == "":
                
            cut_val = parts[-2]

            row = []

            if cut_val in dict.keys():
                row = dict[cut_val]


            f = open(line.strip())
            
            if (db_name == "marker_lib"):
                db_name = db_name + ".18"

            row.append(db_name)

            for ll in f:
                if (ll[0:3] == "tp="):
                    
                    stuff = ll[3:].split()
                        
                    row.append(stuff[0])
                elif (ll[0:3] == "fn="):
                    stuff = ll[3:].split()

                    row.append(stuff[0])


            dict[cut_val] = row
            break
            
        


    
outf = open(out_fn + "all_cnt_out.csv" , "w")
    
vv = dict.keys()
    
vv.sort(key=lambda f: int(f))


for val in vv:

    outf.write(val)
        
    row = dict[val]

    for cell in row:
        
        outf.write("," + cell)
        
    outf.write("\n")

outf.close()

         
    





    
