#!/usr/bin/python

from sys import *

usage = '''
usage: %s rank_fn taxonomy_fn.dat output_fn'

function: outputs a mapping file: tax ID -> rank
''' % argv[0]

if len(argv) != 4 :
  print usage
  exit(1)

out = open(argv[3], 'w')

#build mapping: taxid -> parent
a = open(argv[2])
tid_to_parent = {}
a.readline()
a.readline()
a.readline()
while True :
  x = a.readline()
  y = a.readline()
  if len(x) == 0 or len(y) == 0 :
    break
  t = x.split()
  tid_to_parent[t[0]] = t[-1]

a = open(argv[1])
tid_to_rank = {}
for line in a :
  t = line.split(',')
  t = t[1].split('=')
  taxid = t[1]
  j = line.rfind('\t')
  try :
    assert(j != -1)
    tmp = line[j+1:]
    t = tmp.split(',')
    rank = t[0]
    tid_to_rank[taxid] = rank
  except :
    print line

rank = ''
for x in tid_to_parent.keys() :
 if x != '1' :
  if not tid_to_rank.has_key(x) :
    parent = tid_to_parent[x]
    if not tid_to_rank.has_key(parent) :
      print x,parent
    else :
      rank = tid_to_rank[parent]
  else :
    rank = tid_to_rank[x]

  rank = rank.replace(' ', '_')
  out.write(x + ' ' + rank + '\n')

print len(tid_to_parent)
print len(tid_to_rank)

out.close()
