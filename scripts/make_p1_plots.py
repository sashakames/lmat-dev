import sys

import os


names = ["Wrong Count", "Miss Count", "Unclassified Reads"]
bases = [4, 5, 6]

thresh = ["1", "2", "4", "6", "7"]

for fn in sys.stdin:

    pp = fn.split("/")
    title = pp[-1].strip()


    for nn, bb in zip(names, bases):

        ytxt = ""
        nametxt = ""

        for i in range(4):
            
            if (i == 0 ):
                 y = "y"
                 nametxt = ' name=LO-threshold-' + thresh[1]
            else:
                y = "y" + str(i+1)
                nametxt = nametxt + ' name' + str(i+1) + '=LO-threshold-' + thresh[i+1]
            ytxt = ytxt + " " + y + "=" + str(bb+ ((1+i)*5))

            cmd = '/g/g19/ames4/pl241src/src/pl -prefab lines x=1 ' + ytxt + nametxt + ' xlbl="Min Reads" ylbl="' + nn + '" title="' + title + '" data=' + fn.strip() + ' delim=comma -o ' + fn.strip() + '.' +  str(bb)  + '.eps'

            print cmd
            os.system(cmd)
        
            os.system('ps2pdf -dEPSCrop -dEPSFitPage -dAutoRotatePages=/None ' + fn.strip() + '.' +  str(bb)  + '.eps' )
            os.system('convert -density 100x100 ' + fn.strip() + '.' +  str(bb)  + '.eps '+ fn.strip() + '.' +  str(bb) +   '.png' )
