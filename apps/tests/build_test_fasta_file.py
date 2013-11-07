#!/usr/bin/python

a = open('/usr/mic/post1/metagenomics/ref_sets/microbe2.fasta')

j=-1
while True :
  j+=1
  h = a.readline()
  if len(h) == 0 : break
  if j == 1000 : break;
  i = a.readline()
  assert(len(i) > 0)
  print h,
  print i[:1000]
