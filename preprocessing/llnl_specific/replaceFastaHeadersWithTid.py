#!/usr/bin/python

from sys import *

usage = '''
usage: %s multi output_fn fasta_fn_1 [fasta_fn_2[...]]

function: replaces fasta headers with working tax IDs
''' % (argv[0])

if len(argv) < 4 :
  print usage
  exit(9)


#read tid to header mapping file
a = open(argv[1])
map = {}
for line in a :
  line = line.strip()
  j = line.find('\t')
  assert(j != -1)
  j2 = line.find('\t', j+1)
  assert(j2 != -1)
  j3 = line.find('\t', j2+1)
  assert(j3 != -1)
  j4 = line.find('\t', j3+1)
  assert(j4 != -1)
  header = line[j4+1:]
  working_tid = line[:j]
  map[header] = working_tid


out = open(argv[2], 'w')
is_good = False
seq_count_total = 0
discard = 0
for fn in argv[3:] :
  a = open(fn)
  print 'processing:', fn
  seq_count = 0
  prev_line = ''
  while True :
    if prev_line == '' :
      line = a.readline()
    else :
      line = prev_line
    if len(line) == 0 : break

    if line[0] == '>' :
      line = line.strip()
      if not map.has_key(line) :
        print 'failed to find header:', line
        is_good = False
        discard += 1
        prev_line = a.readline()
        assert(len(prev_line))
      else :
        seq_count += 1
        seq_count_total += 1
        is_good = True
        out.write('>' + map[line] + '\n')
        prev_line = a.readline()
        assert(len(prev_line))
    else :
      if is_good :
        line = line.strip()
        out.write(line)
      prev_line = a.readline()
      if len(prev_line) == 0 : break
      if prev_line[0] == '>' :
        if is_good :
          out.write('\n')
  print '    sequences processed in this file:', seq_count

out.close()
