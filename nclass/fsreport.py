#!/usr/bin/env python

import sys
import pprint

fsfile=sys.argv[1] # fastsummary file
taxtree=sys.argv[2] # ncbi_taxonomy.dat
rankfile=sys.argv[3] # tid to rank mapping
rank_lst_str=sys.argv[4] #ranks of interest
odir=sys.argv[5] #output directory

ranktable={}

a = open(rankfile)
for lin in a :
   vals=lin.split()
   ranktable[vals[0]] = vals[1]

a = open(taxtree)
a.readline()
a.readline()
a.readline()

#read in existing taxonomy
children = {}
names = {}
parent = {}
while True :
  f = a.readline()
  if len(f) == 0 : break
  name = a.readline()
  assert(len(name) > 0)
  t = f.split()
  tid = t[0]
  parent[tid] = t[-1]
  names[tid] = name[:-1]
  children[tid] = {}
  for x in t[2:-1] :
    children[tid][x] = 0

rank_lst=rank_lst_str.split(',')

def getRankTid(rank,tid,ranks,parents) :
   stid=tid
   res=-1
   ## assume this is a plasmid
   if not ranks.has_key(tid) and rank == "plasmid" :
      return tid
   while parent[stid] != stid :
      if  ranks[parent[stid]] == rank :
         res=parent[stid]
         break
      stid = parent[stid]
   return res
         
orig={}
store={}
a = open(fsfile)
for line in a :
   line = line.rstrip()
   t = line.split('\t')
   wrc = t[0]
   count = t[1]
   avg=float(wrc)/float(count)
   taxid = t[2]
   descrip=''
   orig[taxid]=t[3]
   if not parent.has_key(taxid) :
      taxid=-1
      print 'error: failed to find ktaxid', taxid,'for entry:'
      print line
      continue
   else :
      for rank in rank_lst :
         tid=getRankTid(rank,taxid,ranktable,parent)
         if tid == -1 :
            continue
         if not store.has_key(rank) :
            store[rank]={}
         if not store[rank].has_key(tid) :
            store[rank][tid]=[]
         store[rank][tid].append( (taxid,wrc,count) )

for rank in store.keys() :
   fsfileout=odir+fsfile+"."+rank
   print "create fastsummary file ",fsfileout, "for rank=",rank 
   fh=open(fsfileout,"w")
   save=[]
   for tid in store[rank].keys() :
      if orig.has_key(tid) :
         name_str=orig[tid]
      else :
         name_str=names[tid]
      lst=store[rank][tid]
      best_wrc,best_count=-1,-1
      top_strain=-1
      wrc_sum,count_sum=0,0
      for taxid,wrc,count in lst :
         wrc_sum += float(wrc)
         count_sum += int(count)
         if not ranktable.has_key(taxid) :
            ranktable[taxid] = "plasmid"
         if rank == "species" and ranktable[taxid] == "strain" :
            if best_wrc < float(wrc) :
               top_strain=taxid
               best_wrc = float(wrc)
               best_count = count
      strain_info=""
      if top_strain != -1 :
         strain_info = "\t"+str(best_wrc) + "\t"+str(best_count)+"\t"+top_strain+"\t"+orig[top_strain]
      tup=(wrc_sum,count_sum,tid,name_str,strain_info)
      #ouYPt_str=str(wrc_sum)+"\t"+str(count_sum) + "\t" + str(tid)+"\t"+name_str + strain_info
      save.append(tup)
   sval=sorted(save, key=lambda val : val[0],reverse=True)
   for val in sval :
      out_str=str(val[0])+"\t"+str(val[1]) + "\t" + str(val[2])+"\t"+val[3] + val[4]
      fh.write(out_str +"\n")
