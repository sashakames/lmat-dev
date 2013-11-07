import sys

import os


names = ["Wrong Count", "Miss Count", "Unclassified Reads"]
bases = [4, 5, 6]

dbs = ["17", "18", "20"]

fn = sys.argv[1]

pp = fn.split("/")
title = pp[-1].strip()


for nn, bb in zip(names, bases):

    ytxt = ""
    nametxt = ""

    for i in range(3):
            
        if (i == 0 ):
            y = "y"
            nametxt = ' name=microbe2-' + dbs[i]
        else:
            y = "y" + str(i+1)
            nametxt = nametxt + ' name' + str(i+1) + '=microbe2-' + dbs[i]
        ytxt = ytxt + " " + y + "=" + str(bb+ (i*5))

    cmd = '/g/g19/ames4/pl241src/src/pl -prefab lines x=1 ' + ytxt + nametxt + ' xlbl="Min Reads" ylbl="' + nn + '" title="' + title + '" data=' + fn.strip() + ' delim=comma -o ' + fn.strip() + '.' +  str(bb)  + '.eps'

    print cmd
    os.system(cmd)
        
    os.system('ps2pdf -dEPSCrop -dEPSFitPage -dAutoRotatePages=/None ' + fn.strip() + '.' +  str(bb)  + '.eps' )
