#!/usr/bin/python

from sys import *

usage = '''
usage: %s taxonomy_fn 32To16Map_fn output_fn

''' % argv[0]

if len(argv) != 4 :
  print usage
  exit(9)

#read in the mapping
mp = {}
a = open(argv[2])
for line in a :
  t = line.split()
  mp[t[0]] = t[1]
print 'read', len(mp), 'entries from 32->16bit mapping file'

#read in the taxonomy
a = open(argv[1])
a.readline()
a.readline()
a.readline()

parents = {}
names = {}
children = {}

while True :
  line = a.readline()
  if len(line) < 2 : break
  name = a.readline()
  assert(len(name))
  t = line.split()
  tid = t[0]
  parent = t[-1]
  parents[tid] = parent
  names[tid] = name
  children[tid] = {}
  for child in t[2:-1] :
    children[tid][child] = 0

#output the taxonomy subtree, with 16bit ids
out = open(argv[3], 'w')
out.write('#junk\n')
out.write('#junk\n')
out.write('0\n')
ct_nodes = 0
for tid in mp.keys() :
  nid = mp[tid]
  child_ct = 0
  for child in children[tid] :
    if mp.has_key(child) :
      child_ct += 1
  out.write(nid + ' ' + str(child_ct) + ' ')
  for child in children[tid] :
    if mp.has_key(child) :
      out.write(mp[child] + ' ')
  out.write(mp[parents[tid]] + '\n')
  out.write(names[tid])

out.close()
