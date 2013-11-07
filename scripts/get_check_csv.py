import sys

arr = []

dbs = ["marker_lib_17", "marker_lib", "marker_lib_20"]


dict = {}

for d in dbs:
    
#    each entry a row
    dict[d] = {}


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
                
                
                summ_val = parts[-4]
                check_val = parts[-2]
                
                db_dict = dict[db_name]

                row = []

                if check_val in db_dict.keys():
                    row = db_dict[check_val]


                f = open(line.strip())


                for ll in f:
                    if (ll[0:6] == "final:"):
                        
                        llp = ll.strip().split()
                        
                        row.append(summ_val)
                        row.append(llp[2])
                        row.append(llp[4])
                        row.append(llp[6])
                        row.append(llp[8])
                        db_dict[check_val] = row
                        dict[db_name] = db_dict
                break
            
        

for dd in dbs:
    
    outdb = dict[dd]

    if (dd == "marker_lib"):
        dd = dd + ".18"
    
    outf = open(out_fn + dd + ".all_check.csv" , "w")
    
    vv = outdb.keys()
    
    vv.sort(key=lambda f: int(f))


    for val in vv:

        outf.write(val)
        
        row = outdb[val]

        for cell in row:
            
            outf.write("," + cell)
        
        outf.write("\n")

    outf.close()

    print "wrote data for " + dd
        
    





    
