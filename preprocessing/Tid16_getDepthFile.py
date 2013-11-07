#!/usr/bin/python

from sys import *

usage = '''
usage: %s taxonomy_rank_fn 32To16Map_fn output_fn
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

count = 0
out = open(argv[3], 'w')
a = open(argv[1])
for line in a :
  j1 = line.find('taxid=')
  assert(j1 != -1)
  j2 = line.find(',', j1)
  assert(j2 != -1)
  j3 = line.find('ktaxid=')
  assert(j3 != -1)
  j4 = line.find(',', j3)
  assert(j4 != -1)
  taxid = line[j1+6:j2]
  ktaxid = line[j3+7:j4]
  if mp.has_key(taxid) :
    count += 1
    id = mp[taxid]
    replaceme = line[j1:j4]
    r = 'taxid=' + id + ',ktaxid=' + id
    n = line.replace(replaceme, r)
    assert(taxid == ktaxid)
    out.write(n)
out.close()

print 'num IDs in mapping file:', len(mp)
print 'num found in rank file:', count
print '(if not equal, could be an error)'

