#!/usr/bin/python

from sys import *
import os
from common import *
import pprint

Test_fn = 'canonical.fa'

K = 10

os.system('rm -f test.fa.*')

#generate canonical kmers using naive algorithm
kmers = []
a = open(Test_fn)
seq = a.readline()
a.close()
for k in range(0, len(seq)-K) :
    good = True
    for i in range(k,k+K) :
      c = seq[i]
      if not (c == 'a' or c == 'c' or c == 'g' or c == 't') :
        good = False
        break
    if good :
      kmer = seq[k:k+K]
      kmer = canonicalize(kmer)
      kmers.append(kmer)


#generate canoncal kmers using CanoncalEncoder.hpp
cmd = 'genCanonicalKmers ' + Test_fn + ' ' + str(K) + ' test.fa'
runme(cmd)

#read in the kmers
a = open('test.fa')
c_kmers = []
for line in a :
  c_kmers.append(line[:-1])

print len(kmers)
print len(c_kmers)

for j in range(len(kmers)) :
  kmers[j] = kmers[j].upper()

assert(len(kmers) == len(c_kmers))
assert(kmers == c_kmers)
print 'TEST PASSED!'
