
import os
import sys

import numpy



FILES_DIR = sys.argv[1]

NTHREADS_IDX = 1



ofilez = os.popen('ls ' + FILES_DIR +'/out.*.dora5')



arr1 = [4, 8, 12, 16, 20, 24, 32, 40, 64, 80, 120]

arr2 = []

for x in arr1:
   arr2.append([])


for line in ofilez:
       
       tval = 0

       parts = line.strip().split('.')
       f = open(line.rstrip());
       
       
       for ll in f:
          pp = ll.split()
      
          if (pp[1] == "time:"):
             tval = float(pp[-1].rstrip())


       f.close()

       idx = arr1.index(int(parts[NTHREADS_IDX]))
          
       if (idx > -1):
          arr2[idx].append(tval)
   


for x, y in zip(arr1, arr2):


    print str(x) +  " " + str(numpy.mean(y)) + " " + str(numpy.std( y ))


