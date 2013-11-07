#!/usr/bin/python

import sys
import os

usage = '''
usage: %s input_fn output_dir <int>
where: <int> is the approx. size, in G, of the partitioned files.
function: partitions a fasta file into smaller pieces.
Note: output directory will be created if it doesn't exist.
      .int files have a \\n between header and seq so use on org .fasts
''' % sys.argv[0]

if len(sys.argv) != 4 :
  print usage
  sys.exit(-1)


a = open(sys.argv[1])
os.system('mkdir -p '+ sys.argv[2])

output_file_size = float(sys.argv[3]) *1000000000

file_ext = 0
nm = sys.argv[1]
j = sys.argv[1].rfind('/')
if j != -1 : nm = sys.argv[1][j+1:]
out = open(sys.argv[2] + '/' + nm + '.part.' + str(file_ext), 'w')
ct =  0
while True :
  if ct >= output_file_size :
    out.close()
    file_ext += 1
    out = open(sys.argv[2] + '/' + nm + '.part.' + str(file_ext), 'w')

    ct = 0
  header = a.readline()
  seq = a.readline()
  if len(header) == 0 :
    assert(len(seq) == 0)
    break

  #this will check that sequence is contained on a single line
  assert(header[0] == '>')

  out.write(header)
  out.write(seq)
  ct += len(seq)
  ct += len(header)
out.close()



