import sys

import os


names = ["True Positive Pct", "False Negative Pct", ]
bases = [3, 4]

dbs = ["marker_lib_17", "marker_lib.18", "marker_lib.20"]

fn = sys.argv[1]

pp = fn.split("/")
title = pp[-1].strip()


for nn, bb in zip(names, bases):

    ytxt = ""
    nametxt = ""

    for i in range(3):
            
        if (i == 0 ):
            y = "y"
            nametxt = " name=" +dbs[i]
        else:
            y = "y" + str(i+1)
            nametxt = nametxt + ' name' + str(i+1) + '=' + dbs[i]
        ytxt = ytxt + " " + y + "=" + str(bb+ (i*3))

    cmd = '/g/g19/ames4/pl241src/src/pl -prefab lines x=1 ' + ytxt + nametxt + ' xlbl="LO Thresh." ylbl="' + nn + '" title="' + title + '" data=' + fn.strip() + ' delim=comma -o ' + fn.strip() + '.' +  str(bb)  + '.eps'

    print cmd
    os.system(cmd)
        
#    os.system('ps2pdf -dEPSCrop -dEPSFitPage -dAutoRotatePages=/None ' + fn.strip() + '.' +  str(bb)  + '.eps' )
#    os.system('convert -density 100x100 ' + fn.strip() + '.' +  str(bb)  + '.eps '+ fn.strip() + '.' +  str(bb) +   '.png' )
