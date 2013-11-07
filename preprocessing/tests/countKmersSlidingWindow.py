#!/usr/bin/python

from sys import *
import os
from common import *
import pprint

K = 10

#generate canonical kmers using naive algorithm
kmers = {}
a = open(argv[1])
ct = 0
while True :
  b = a.readline()
  if len(b) == 0 : break;
  seq = a.readline()
  assert(len(seq) > 0)
  for k in range(0, len(seq)-K) :
    good = True
    for i in range(k,k+K) :
      c = seq[i].lower()
      if not (c == 'a' or c == 'c' or c == 'g' or c == 't') :
        good = False
        break
    if good :
      kmer = seq[k:k+K]
      kmer = canonicalize(kmer)
      kmers[kmer] = 0
      ct += 1
      if ct % 1000000 == 0 : print ct/1000000, 'M'

print 'kmer count:', len(kmers)
