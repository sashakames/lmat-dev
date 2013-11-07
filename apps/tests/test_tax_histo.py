#!/usr/bin/python

import sys
import os

'''
runs build_header_lookup_table.py on vbonly_100;
runs jellylist on test.fa.int;
runs tax_histo_fast_limited and tax_histo_fast_limited_old on output;
tests that output from  tax_histo_fast_limited_old contains duplicate tax IDs;
tests that output from  tax_histo_fast_limited does not contain duplicate tax IDs;
tests that output from  tax_histo_fast_limited_old and tax_histo_limited have the same
set of tax IDs for each kmer

'''

cur = os.getcwd()

def runme(cmd) :
  r = os.system(cmd)
  if r != 0 :
    sys.stderr.write('TEST FAILED\n')
    sys.exit(-1)
  else :
    sys.stderr.write('PASSED\n')


cmd = 'rm -rf tmp; mkdir tmp'
sys.stderr.write('erasing and recreating tmp directory\n')
runme(cmd)
cmd = 'rm -rf tmp; mkdir tmp'
runme(cmd)

cmd = 'cp test.fa tmp'
sys.stderr.write('copying test.fa to tmp\n')
runme(cmd)

cmd = '../../preprocessing/build_header_lookup_table.py test.fa tmp/test.fa_table'
sys.stderr.write('running: ' + cmd + '\n')
runme(cmd)

cmd = 'mv test.fa.int tmp'
sys.stderr.write('mv test.fa.int tmp\n\n')
runme(cmd)

cmd = '../../preprocessing/jellylist_bin -h 1000000 -t 32 -k 18  -loc 0  tmp/test.fa.int'
sys.stderr.write('running jellylist: ' + cmd + '\n\n')
runme(cmd)

cmd = "../tax_histo_fast_limited -o tmp/output_new -d tmp/test.fa.int.bin.no_loc.18 -a tmp/test.fa_table -g ../../runtime_inputs/gid_to_tid_all.txt -p 3 -t ../../runtime_inputs/kpath_taxonomy.dat"
sys.stderr.write('running tax_histo_fast_limited\n')
runme(cmd)

cmd = "../tax_histo_fast_limited_old -o tmp/output_old -d tmp/test.fa.int.bin.no_loc.18 -a tmp/test.fa_table -g ../../runtime_inputs/gid_to_tid_all.txt -p 3 -t ../../runtime_inputs/kpath_taxonomy.dat"
sys.stderr.write('running tax_histo_fast_limited_old\n')
runme(cmd)

cmd = "testTaxHisto tmp/output_old tmp/output_new"
sys.stderr.write('testing that tax IDs are identical sets for tax_histo_new and tax_histo_old')
runme(cmd)
