#!/usr/bin/python

from sys import *
from common import *
import os

if len(argv)!=2 :
  print 'usage: %s output_fn_from_tax_histo' % argv[0]
  exit(1)

cmd = '../convertTaxHisto_ascii2bin ' + argv[1] + ' .'
runme(cmd)
name = argv[1]
z = name.rfind('/')
if z != -1 :
  name = name[z+1:]


cmd = '../convertTaxHisto_bin2ascii ' + name + '.bin output_ascii2bin2ascii'
runme(cmd)

print 'reading:', argv[1]
a = open(argv[1]).readlines()
print 'reading: output_ascii2bin2ascii'
b = open('output_ascii2bin2ascii').readlines()
print len(a), len(b)
#assert(len(a)-2 == len(b))
for j in range(2, len(a)) :
  t = a[j].split()
  t2 = b[j-2].split()
  x1 = {}
  x2 = {}
  for z in t : x1[z] = 0
  for z in t2 : x2[z] = 0
  if x1 != x1 :
    print 'failed for j=', j
    exit(1)

print 'PASSED!'
