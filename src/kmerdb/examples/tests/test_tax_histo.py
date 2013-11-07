#!/usr/bin/python

import sys

if len(sys.argv) != 3 :
  print 'usage: test_tax_histo.py input1 input2'
  print 'where: input1,2 are the output files from two of the tax_histo programs.'
  sys.exit(-1)

print 'test 1: ======================================================='
a = open(sys.argv[1]).readlines()
b = open(sys.argv[2]).readlines()

s = {}
for j in range(2, len(a)) :
  t = a[j].split()
  kmer = t[0]
  for x in range(3, len(t), 2) :
    try :
      s[(t[x],t[x+1], kmer)] = 0
    except :  
      print 'error on line', j

nf1 = 0
for j in range(2, len(b)) :
  t = b[j].split()
  kmer = t[0]
  for x in range(3, len(t), 2) :
    if not s.has_key((t[x],t[x+1], kmer)) :
      print '1. not found:', t[x],t[x+1],kmer,'on line',j
      nf1 += 1

print 'test 2: ======================================================='

s = {}
for j in range(2, len(b)) :
  t = b[j].split()
  kmer = t[0]
  for x in range(3, len(t), 2) :
    try :
      s[(t[x],t[x+1], kmer)] = 0
    except :  
      print 'error on line', j

nf2 = 0
for j in range(2, len(a)) :
  t = a[j].split()
  kmer = t[0]
  for x in range(3, len(t), 2) :
    if not s.has_key((t[x],t[x+1], kmer)) :
      print '1. not found:', t[x],t[x+1],kmer,'on line',j
      nf1 += 1

print
print 'not found in test 1:', nf1
print 'not found in test 2:', nf2
