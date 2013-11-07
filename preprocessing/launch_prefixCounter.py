#!/usr/bin/env python

from sys import *
import os

header = '''
#MSUB -S /bin/tcsh
#MSUB -l partition=sierra
#MSUB -A prophecy
#MSUB -q pbatch
#MSUB -l nodes=1
#MSUB -l walltime=30:00
#

'''

if len(argv) != 5 :
  print 'usage:', argv[0], ' prefix_size input_fasta_file (preprocessed) output_basename ksize'
  print 'example: ', argv[0], ' 4 Bacteria.genes.fa.int Bacteria.kmerDB 20'
  print 'WARNING, SCRIPTS HARDCODED FOR MACHINE/BANK/QUEUE (can change this eventually as needed)
  exit(1)

psize = int(argv[1])
fasta_in = argv[2]
out_name = argv[3]
ksize = int(argv[4])
num_jobs=pow(4,psize)
prog='/p/lscratchc/allen99/metag/metag_repo/rel/LMAT-1.1/src/kmerPrefixCounter'
for j in range(0,num_jobs) :
  out = open('msub.' + str(j), 'w')
  out.write(header)
  out.write(prog + ' -k ' + str(ksize) + ' -i ' + fasta_in + ' -o ' + out_name + ' -l ' + str(psize) + ' -f ' + str(j) + '\n')
  out.close()

