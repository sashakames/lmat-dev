#!/usr/bin/python

from sys import *

usage = '''
usage: %s <read_label_output_fn> 32To16Map_fn output_fn

''' % argv[0]

if len(argv) != 4 :
  print usage
  exit(9)

#read in the mapping
mp = {}
a = open(argv[2])
for line in a :
  t = line.split()
  mp[t[0]] = t[1]
print 'read', len(mp), 'entries from 32->16bit mapping file'
a.close()

a = open(argv[1])
out = open(argv[3], 'w')

ln = 0
for line in a :
  ln += 1
  t = line.split('\t')
  t2 = []
  t3 = t[3].split()
  for j in range(0, len(t3), 2) :
    if not mp.has_key(t3[j]) :
      print 'ERROR; tax ID', t3[j], 'was not found in mapping file; error on line', ln, 'of read_label file, for j=', j
      exit(9)
    else :
      t2.append(mp[t3[j]])
      t2.append(t3[j+1])

  t4 = t[4].split()
  t44 = []
  if not mp.has_key(t4[0]) :
      print 'ERROR; tax ID', t4[0], 'was not found in mapping file; error on line', ln, 'of read_label file, for final field'
      exit(9)
  else :
      t44.append(mp[t4[0]])
      t44.append(t4[1])
      t44.append(t4[2])

  tt = ' '.join(t2)
  ttt = ' '.join(t44)
  out.write(t[0] + '\t' + t[1] + '\t' + t[2] + '\t' + tt + '\t' + ttt + '\n')

a.close()
out.close()

'''
> r1.1 
> |SOURCES={GI=384055883,fw,158600-158700}|ERRORS={80:G,86:C}|SOURCE_1="Acet
> obacter pasteurianus IFO 3283-22"
> (e0ef0893d4de1486c414f36d1ad36b23fd275fa7) \t
> CGCCCAGCGGATGAAAAGCTGCGCAGCCCGCCACCACGGACCCAGTTCAGCGGGGATACTCCGCCATTTCGCGC
> CATTCTGATGACCCCAGAAAATCGCT \t 0.0421036 1.9824 81 \t 940284 3.92607
> 10002616 3.83471 634455 3.83471 634452 3.83471 1004836 3.83471 33995
> 3.75869 89583 3.75869 433 3.75869 10001041 3.63941 10000423 3.63941
> 10000422 3.63941 10000775 3.63941 10001901 3.63941 10000690 3.63941
> 10000308 3.63941 634457 3.63941 634458 3.63941 940265 3.63941 945681
> 3.63941 634177 3.63941 438 3.63941 434 3.63941 204441 3.61169 634453
> 3.47203 634456 3.47203 634459 3.47203 481145 3.47203 65959 3.47203 28448
> 3.47203 634454 3.19357 28211 2.88274 887700 2.74426 120045 2.74426 542
> 2.74426 541 2.74426 627344 2.59596 435 2.54896 151157 2.54896 1224 2.5164
> 1006554 2.45016 41297 2.10312 204457 2.10312 940285 1.94032 940286
> 1.91309 131567 1.7385 1 1.7385 2 1.72518 65958 1.7178 290633 1.60443 442
> 1.60443 441 1.40913 10001981 1.31777 10002449 1.20618 622759 0.97037
> 940282 0.775071 940283 0.396794 264203 0.361415 120044 0.234088 \t 28211
> 3.92607 MultiMatch

'''
