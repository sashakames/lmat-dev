#!/usr/bin/python

from sys import *
import os
from common import *
import pprint


usage = '''
usage: %s mer_len fasta_fn output_from_prefixSorter
''' % argv[0]

if len(argv) != 4 :
  print usage
  exit(-1)

#read in output from kmerPrefixCounterSort
a = open(argv[3])
s = {}
for line in a :
  t = line.split()
  assert(not s.has_key(t[0]))
  s[t[0]] = {}
  for id in t[1:] :
    s[t[0]][id] = 0

#generate canonical kmers using naive algorithm
K = int(argv[1])
kmers = {}
a = open(argv[2])
tid = -1
while True :
  seq = a.readline()
  if len(seq) == 0 : break;
  if seq[0] == '>' :
    tid = seq[1:-1]
    print 'tid:', tid
  else :
    #print 'len:',len(seq)
    #print seq.lower(),
    for k in range(0, len(seq)-K) :
      good = True
      for i in range(k,k+K) :
        c = seq[i].lower()
        if not (c == 'a' or c == 'c' or c == 'g' or c == 't') :
          good = False
          break
      if good :
        kmer = seq[k:k+K]
        #print k,kmer.lower(),
        kmer2 = canonicalize(kmer)
        #print kmer2
        if not kmers.has_key(kmer2) : kmers[kmer2] = {}
        kmers[kmer2][tid] = 0

#test if kmer sets are equal
bad = 0
for kmer in s.keys() :
  if not kmers.has_key(kmer) :
    bad += 1
print 'in s, not in kmers:', bad

bad = 0
for kmer in kmers.keys() :
  if not s.has_key(kmer) :
    bad += 1
print 'in kmers, not in s:', bad

k1 = s.keys()
k2 = kmers.keys()
k1.sort()
k2.sort()
#for j in range(20) : print k1[j], k2[j], k1[j] == k2[j]
if k1 == k2 : print 'kmer sets are equal'
else :
  print 'kmer sets are NOT equal; FAILED!'
  exit(-1)



#test if tax ID sets are equal
ct = 0
for kmer in k1 :
  ct += 1
  if not s[kmer] == kmers[kmer] :
    print 'error on kmer number', ct
    print 'kmer:', kmer
    print 'tax IDs, 1:', s[kmer]
    print 'tax IDs, 2:', kmers[kmer]
    exit(9)
print 'tax ID sets are equal'
