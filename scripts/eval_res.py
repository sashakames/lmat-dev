import sys

import os

#databases = ["kML.18.noOvr.no_prune.16bit.db", "kML.18.full.16bit.db", "kML.20mer.16bit.g1000.db"  ]
databases = ["kML.18.noOvr.no_prune.16bit.db",  "kML.20mer.16bit.g1000.db"  ]
#[ "kML.18.full.16bit.db", 
#sizes = [1000, 9000, 27000, 81000, 2187000]
sizes = [2000000]
positions = ["head", "tail"]
thresholds = [20, 40 ]


#print "Input\tDatabase\tThreshold\tSize\tPosition\tTotal\tStrain match\tstrainotherspecies\tstrainothergenus\tWrong strain\tSpecies match\tOther species\tWrong species\tGenus match\tWrong genus"

bigarr = []

for line in sys.stdin:

    if line[0:3] == "fc.":
        continue

    name = line.strip()

    name_parts = name.split("_")


    strain = name.replace("_", " ")
    species = name_parts[0] + " " + name_parts[1]
    genus = name_parts[0]

    for db in databases:

        for thresh in thresholds:
            
            for sz in sizes:
            
                for pos in positions:



                    fname = sys.argv[1] + "/dir.", thresh, ".", name,".", db,".", sz, "/", name, ".", sz, ".", pos, ".", db, ".", thresh, ".lo.rl_output.0.5.fastsummary"

                    fname = ''.join(map(str, fname))
                    #                    print fname
                    if not os.path.exists(fname):
                        continue

                    # initialize counters

                    f = open(fname)
                    
                    total = 0
                    strainmatch = 0
                    strainotherspecies = 0
                    strainothergenus = 0
                    strainwrong = 0

                    speciesmatch = 0

                    speciesother = 0
                    specieswrong = 0
                    
                    genusmatch = 0
                    genuswrong = 0


                    for result in f:
                        

                        parts = result.split("\t")
                        count = int(parts[1])
                        
                        total += count

                        if parts[3] == "no rank,cellular organisms\n":
                            continue

                        pair = parts[3].replace("[","").replace("]","").split(",")
                        
                        rank = pair[0]

                        namewords = pair[1].strip().split(" ")

                        # check for strains
                        if (rank == "no rank" and len(namewords) > 1):
                            
                            if (strain == pair[1].strip()):
                                strainmatch += count
                            elif (name_parts[0] == namewords[0]):

                                if (name_parts[1] == namewords[1]):
                                    
                                    strainotherspecies += count
                                else:
                                    strainothergenus += count
                            else:
                                strainwrong += count

                        elif (rank == "species"):
                            
                            if (name_parts[0] == namewords[0]):

                                if (name_parts[1] == namewords[1]):
                                    
                                    speciesmatch += count
                                else:
                                    speciesother += count
                            else:
                                specieswrong += count
                        elif (rank == "species group"):
                            if (genus == namewords[0]):
                                genusmatch += count;
                            else:
                                genuswrong += count;
                        elif (rank == "genus"):
                            if (genus == pair[1].strip()):
                                genusmatch += count
                            else:
                                genuswrong += count

                    print name, "\t", db, "\t", thresh, "\t", sz, "\t", pos, "\t", total, "\t", strainmatch, "\t", strainotherspecies, "\t", strainothergenus, "\t", strainwrong, "\t", speciesmatch, "\t", speciesother, "\t", specieswrong, "\t", genusmatch, "\t", genuswrong


#                    bigarr.append([ name, db, thresh, sz, pos, total, strainmatch, strainotherspecies, strainothergenus, strainwrong, speciesmatch, speciesother, specieswrong, genusmatch, genuswrong])






    
