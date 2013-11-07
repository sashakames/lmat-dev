import random
import os
import sys

def runme(cmd) :
  sys.stderr.write('running: ' + cmd + '\n')
  r = os.system(cmd)
  if r != 0 :
    sys.stderr.write('TEST FAILED\n')
    sys.exit(-1)

def canonicalize(kmer) :
  kmer2 = kmer.lower()
  rev = rc(kmer)
  if kmer2 < rev : return kmer2
  return rev


def rc(kmer) :
  t = []
  for c in kmer :
    c = c.lower() 
    if c == 'a' : t.append('t')
    elif c == 't' : t.append('a')
    elif c == 'c' : t.append('g')
    elif c == 'g' : t.append('c')
    else :
      print 'error; not a, c, g, or t:', c
      assert(False)
  km = ''.join(t)
  return km[::-1]

def generateSeqs(number, length) :
  out = []
  j = 0
  for x in range(number) :
    header = '>' + str(j);
    j += 1
    out.append(header)
    z = []
    for y in range(length) :
      s = random.randint(0,3)
      if s == 0 : z.append('a')
      elif s == 1 : z.append('c')
      elif s == 2 : z.append('g')
      elif s == 3 : z.append('t')
      else :
        assert(False)
    out.append(''.join(z))
  return out

def writeSeqs(seqs, num, last, fn) :
  out = open(fn, 'w')
  for j in range(first*2, (first+num)*2) :
    out.write(seqs[j] + '\n');
  out.close()
