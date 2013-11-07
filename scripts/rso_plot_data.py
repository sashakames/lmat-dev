
import os
import sys

import numpy



FILES_DIR = sys.argv[1]

NTHREADS_IDX = 1

MACH_NAME = sys.argv[2]






#arr1 = [4, 8, 12, 16, 20, 24, 32, 40, 64, 80, 120]
#arr1 = [.1, .5, 1, 2, 5, 10, 20]
arr2 = []
#arr1 =[.01, .05, .1, .2, 15, 20, 30, 40]


argsfile = open(FILES_DIR + '/argsfile.' + MACH_NAME)

line = argsfile.readline()

line = argsfile.readline()

parts = line.strip().split()

arr1 = []

for x in parts[1:]:

   arr1.append(float(x))

#arr1 = [0, .01, .05, .1, .2, .5, 1, 2, 5, 10, 15, 20, 30, 40]

argsfile.close()

for x in arr1:
   arr2.append([])


ofilez = os.popen('ls ' + FILES_DIR +'/out.*.' + MACH_NAME)




for line in ofilez:
       
       tval = 0

       ss = line.replace("..", ".;").replace(".",":").replace(';','.')
       
       parts = ss.strip().split(':')
       f = open(line.rstrip());
       
       
       for ll in f:
          

          pp = ll.split()
      
          if (pp[1] == "time:"):
             tval = float(pp[-1].rstrip())


       f.close()

       idx = arr1.index(float(parts[NTHREADS_IDX]))
          
       if (idx > -1):
          arr2[idx].append(tval)
   


for x, y in zip(arr1, arr2):


    print str(x) +  " " + str(numpy.mean(y)) + " " + str(numpy.std( y ))


