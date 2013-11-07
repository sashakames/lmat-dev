#!/usr/bin/python

from sys import *
import os
from common import *
import pprint

Test_fn = 'test.fa'
Test_fn = 'ten_seqs.fa'

K = 10

os.system('rm -f test.fa.*')

kmers = {}
a = open(Test_fn).readlines()
for j in range(0,len(a),2) :
  gid = a[j][1:-1]
  for k in range(0, len(a[j+1])-K) :
    good = True
    for i in range(k,k+K) :
      c = a[j+1][i]
      if not (c == 'a' or c == 'c' or c == 'g' or c == 't') :
        good = False
        break
    if good :
      kmer = a[j+1][k:k+K]
      kmer = canonicalize(kmer)
      #print kmer, kmer2
      if not kmers.has_key(kmer) :
        kmers[kmer] = {}
      kmers[kmer][gid] = 0


'''
code to ensure the test fails:
1. add a kmer to the first set that isn't in the second set
kmers['atatatatat'] = 0

2. add a gid to a kmer in the first set, which isn't in the second set
for k in kmers :
  kmers[k][11] = 0
  break
'''

      
#pprint.pprint(kmers)
print 'kmer count:', len(kmers)

cmd = 'mpirun -np 4 ../kmerPrefixCounter_MPI -i ' + Test_fn + ' -k ' + str(K) + ' -o test.fa -l 1'
print 'running:', cmd
os.system(cmd)

kmers2 = {}
tt = 0;
for j in range(4) :
  cmd = 'dumpBinaryKmers test.fa.' + str(j) + ' > tmp'
  os.system(cmd)
  a = open('tmp')
  ct = 0;
  for line in a :
    ct += 1
    tt += 1
    t = line.split()
    kmer = t[1]
    if kmer == 'ccggccgcgc' : print 'yyy',line, len(t)
    if not kmers2.has_key(kmer) :
      kmers2[kmer] = {}
    for gid in t[2:] :
      kmers2[kmer][gid] = 0  
  print 'read', ct, 'kmers from','dumpBinaryKmers test.fa.' + str(j)
print 'total kmers from file:', tt

b = open('out.kmers','w')
for km in kmers.keys() : b.write(km+'\n')
b.close()
b = open('out.kmers2','w')
for km in kmers2.keys() : b.write(km+'\n')
b.close()
  

print
print
print 'len(kmers.keys():', len(kmers.keys())
print 'len(kmers2.keys():', len(kmers2.keys())

print 
c = 0
print 'in kmers2, not in kmers:'
for km in kmers2.keys() :
  if not kmers.has_key(km) :
    print 'missing:', km
    c += 1
print
if c == 0 : print 'PASSED'
else : print 'FAILED'

c = 0
print 'in kmers, not in kmers2:'
for km in kmers.keys() :
  if not kmers2.has_key(km) :
    k2 = rc(km)
    if not kmers2.has_key(km) :
      print 'missing:', km, k2
      c += 1
print
if c == 0 : print 'PASSED'
else : print 'FAILED'
print

print 'testing that keys match (duplicates the above)'
assert(len(kmers.keys()) == len(kmers2.keys()))

print 'testing that gids match'
failed = False
for kmer in kmers.keys() :
  if not kmers2.has_key(kmer) :
    print 'in kmers, not in kmers2:', kmer
    failed = True
    break
  if kmers[kmer] != kmers2[kmer] :
    print 'error: kmers[', kmer,']:', kmers[kmer], 'kmers2[kmer]',kmers2[kmer]
if failed : print 'FAILED'
else : print 'PASSED'

ca = 0
cc = 0
cg = 0
ct = 0
for k in kmers.keys() :
  if k[0] == 'a' : ca += 1
  elif k[0] == 'c' : cc += 1
  elif k[0] == 'g' : cg += 1
  elif k[0] == 't' : ct += 1
print 
print 'a',ca
print 'c',cc
print 'g',cg
print 't',ct
