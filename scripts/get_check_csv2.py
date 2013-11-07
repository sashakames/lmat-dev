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

    if (parts[-1] == "check" and parts[-3] == "summary"):

        db_name = ""

        for pp in parts:

            for dd in dbs:

                if pp == dd:
                    
                    db_name = pp
            
            if not db_name == "":
                
                

                check_val = parts[-2]
                


                row = []

                if check_val in dict.keys():
                    row = dict[check_val]


                f = open(line.strip())


                for ll in f:
                    if (ll[0:6] == "final:"):
                        
                        llp = ll.strip().split()
                        if (db_name == "marker_lib"):
                            db_name = db_name + ".18"
                        
                        row.append(db_name)
                        row.append(llp[2])
                        row.append(llp[4])
                        row.append(llp[6])
                        row.append(llp[8])
                        dict[check_val] = row

                break
            
        



    
outf = open(out_fn + sys.argv[2]  ".all_check.csv" , "w")
    
vv = dict.keys()
    
vv.sort(key=lambda f: int(f))


for val in vv:

    outf.write(val)
    
    row = dict[val]
        
    for cell in row:
            
        outf.write("," + cell)
            
        outf.write("\n")

outf.close()

        
    





    
