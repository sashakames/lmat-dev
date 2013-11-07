#!/usr/bin/python
from sys import *

usage = '''
usage: %s input_taxonomy_file input_32_16_mapping output_fn

example:
build16bitTaxonomy.py ../runtime_inputs/ncbi_taxonomy.segment_20130204.dat ../runtime_inputs/ncbi_tid_32_to_16_map_20130204.txt ../runtime_inputs/ncbi_taxonomy.segment_20130204_16bit.dat
''' % argv[0]

if len(argv) != 4 :
  print usage
  exit(1)

print 'opening mapping file:', argv[2]
mp = {}
a = open(argv[2])
for line in a :
  t = line.split()
  mp[t[0]] = t[1]


print 'opening output file:', argv[3]
out = open(argv[3], 'w')
out.write('#junk\n')
out.write('#junk\n')
out.write('0\n')

print 'opening taxonomy file:', argv[1]
a = open(argv[1]).readlines()
for j in range(3, len(a), 2) :
  t = a[j].split()
  name = a[j+1]
  parent = t[-1]
  tid = t[0]
  if mp.has_key(t[0]) :
    if not mp.has_key(parent) :
      print 'no parent for tid:', tid
    else :
      parent = mp[parent]
      c = {}
      assert(len(t[2:-1]) == int(t[1]))
      for child in t[2:-1] :
        if mp.has_key(child) :
          c[child] = 0
      out.write(mp[tid] + ' ' + str(len(c)) + ' ');    
      for child in c :
        out.write(mp[child] + ' ')
      out.write(parent +'\n')  
      out.write(name)

out.close()
