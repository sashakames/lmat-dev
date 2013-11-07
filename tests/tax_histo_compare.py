#!/usr/bin/python

import sys

a = open(sys.argv[1]).readlines()
b = open(sys.argv[2]).readlines()
assert(len(a) == len(b))

for j in range(len(a)) :
  t1 = a[j].split()
  t2 = b[j].split()
  s1 = {}
  assert(len(t1) == len(t2))

  if len(t1) == 3:
    assert(t1[2] == t2[2])
  else :
    for z in range(2, len(t1)) :
      y1 = t1[z].split(',')
      y2 = t2[z].split(',')
      assert(y1[0] == y2[0] and y1[1] == y2[1])
    

