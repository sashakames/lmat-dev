#!/usr/bin/python

import sys

usage = '''
usage: %s gid_to_tid_all_input table_fn_input  output_dir gid_to_tid_output not_found_fn_output
''' % sys.argv[0]

if len(sys.argv) != 6 :
  print usage
  for j in range(len(sys.argv)) : print j, sys.argv[j]
  sys.exit(-1)

g2t = {}
a = open(sys.argv[1])
for line in a :
  t = line.split()
  g2t[t[0]] = t[1]

a = open(sys.argv[2]).readlines()
out = open(sys.argv[3]+'/' + sys.argv[4], 'w')
nf = open(sys.argv[3] + '/' + sys.argv[5], 'w')
nf.write('genome IDs that failed to map to tax node IDs:\n\n')

for j in range(0, len(a), 2) :
  gid = a[j][:-1]
  if g2t.has_key(gid) :
    out.write('%s %s\n' % (gid,g2t[gid]))
  else :
    nf.write('  %s ::%s\n' % (gid,a[j+1]))

out.close()
nf.close()
