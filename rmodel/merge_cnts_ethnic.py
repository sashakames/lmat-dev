#!/usr/bin/env python

from sys import *
import pdb

usage = '''
usage: %s threshold lmat_lineage_file taxtree output_strain_select_version output_near_neighbor_map remapped ids ktax2ncbitax (optional)
where: input_taxonomy is in nodes.dmp format
''' % argv[0]

if len(argv) != 8 :
  print "args: ",len(argv)
  print usage
  exit(1)

mag_diff=100
tax_obs=argv[1]
ncbi_taxonomy=argv[2]
ncbi_rank =argv[3]
min_obs=int(argv[4])
print "min_obs:",min_obs
thc_file=argv[5]
output_file=argv[6]
num_bins=int(argv[7])
ranks = {}
tax_hist_cnt={}
thc_strm = open(thc_file)
for line in thc_strm :
   line=line.rstrip()
   vals=line.split()
   tid=int(vals[0])
   cnt=int(vals[1])
   tax_hist_cnt.setdefault(tid,cnt)
    
a = open(ncbi_rank)
for line in a :
  line=line.rstrip()
  tid,rank=line.split()
  ranks[int(tid)] = rank

ranks.setdefault(1,"life")

#get map: tax_id -> parent
a = open(ncbi_taxonomy)
names={}
parents = {}
read_cnt={}
children={}
while True :
  line1=a.readline()
  if not line1 :
     break
  if line1[0] == '#' :
     a.readline()
     a.readline()
     continue
  vals=line1.split()
  tid=int(vals[0])
  num_child=int(vals[1])
  name=a.readline() # ignore
  #names[tid] = name.rstrip()
  parent = int(vals[len(vals)-1].rstrip())
  parents[tid] = parent
  children[tid]={}
  for it in range(0,num_child,1) :
      assert it+2 < len(vals)
      child=vals[it+2] 
      children[tid][child] = 0

obs_val={}
tax_obs_strm = open(tax_obs)
for val in tax_obs_strm :
   val=val.rstrip()
   t=val.split()
   id=int(t.pop(0))
   obs_val.setdefault(id,t)
     
cnt1=0

print "how much",cnt1
cnt=0
out_fh=open(output_file,"w")
out_fh.write(str(num_bins)+"\n")
qlst =  tax_hist_cnt.keys()
print qlst
once={}
save={}
saveid=[]
warn_cnt=0
for tid in qlst :
   if once.has_key(tid) :
      continue
   once.setdefault(tid,1)
   if not parents.has_key(tid) :
      continue
   curr_tid = tid
   use_val = obs_val[tid]
   tid_kcnt=0
   if tax_hist_cnt.has_key(tid) :
      tid_kcnt=tax_hist_cnt[tid]

   isOther=False
   verbose=False
   (rval_pcnt,rval_kcnt,rval_obs)=([0]*num_bins,[0]*num_bins,[0]*num_bins)
   (rval_pcnt1,rval_kcnt1,rval_obs1)=([1.0]*num_bins,[0]*num_bins,[0]*num_bins)
   close_match=[-1]*num_bins
   fndMatch=False
   cnt=0
   for it in range(0,len(use_val),2) :
      rval_pcnt[cnt] = float(use_val[it])
      rval_kcnt[cnt] = tid_kcnt
      rval_obs[cnt] = float(use_val[it+1])
      cnt+=1

   use_rank = ranks[curr_tid]
   str_out=str(tid)+" "+str(use_rank)+"-"+str(curr_tid)
   save_rit,save_fit=-1,-1
   not_min=False
   for it in range(len(rval_pcnt)) :
      if int(rval_obs[it]) < min_obs : 
         not_min=True
         for rit in range(it-1,-1,-1):
            if int(rval_obs[rit]) >= min_obs :
               save_rit=rit
               break
         for fit in range(it+1,len(rval_pcnt),1):
            if int(rval_obs[fit]) >= min_obs :
               save_fit=fit
               break
         if save_rit >= 0 :
            d1=abs(it-save_rit)
         else :
            d1=num_bins+1
         if save_fit >= 0 :
            d2=abs(it-save_fit)
         else :
            d2=num_bins+1
         if d1 <= d2 and save_rit != -1 :
            rval_pcnt[it] = rval_pcnt[save_rit]
         elif save_fit != -1 :
            rval_pcnt[it] = rval_pcnt[save_fit]

   for it in range(len(rval_pcnt)) :
      str_out += " " + str(rval_obs[it])+" "+str(rval_pcnt[it])+" " +str(rval_kcnt[it])
   if save_rit == -1 and save_fit == -1 and not_min :
      warn_cnt+=1
   out_fh.write(str_out+"\n")
   cnt+=1
