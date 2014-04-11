#!/usr/bin/env python

import sys,os
import pprint

fsfile=sys.argv[1] # fastsummary file
taxtree=sys.argv[2] # ncbi_taxonomy.dat
rankfile=sys.argv[3] # tid to rank mapping
rank_lst_str=sys.argv[4] #ranks of interest
odir=sys.argv[5] #output directory
gsfile=""
if len(sys.argv) >= 6 :
   gsfile=sys.argv[6]

#debug_str="[p1:="+fsfile+"]\n[p2:"+taxtree+"]\n[p3:"+rankfile+"]\n[p4:"+rank_lst_str+"]\n[p5:"+odir+"]"
#print "debug",debug_str
ranktable={}

print "open1",rankfile
a = open(rankfile)
for lin in a :
   vals=lin.split()
   ranktable[vals[0]] = vals[1]

print "open2",taxtree
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
      #if  ranks[parent[stid]] == rank :
      if  ranks.has_key(stid) and ranks[stid] == rank :
         res=stid
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
      if taxid != 1 :
         taxid=-1
         print 'warning: did not find parent id for node (ignore)', taxid,'for entry:'
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
gene_store={}
if gsfile != "" : 
   a = open(gsfile)
   for line in a :
      line = line.rstrip()
      t = line.split('\t')
      rc = t[0]
      taxid = t[1]
      type=t[6]
      ### track how many reads were mapped to rRNA
      if type != "rRNA" :
         continue
      
      if not parent.has_key(taxid) :
         taxid=-1
         print 'warning: did not parent node for', taxid,' entry: (rRNA gene not counted)'
         print line
         continue
      else :
         for rank in rank_lst :
            tid=getRankTid(rank,taxid,ranktable,parent)
            if tid == -1 :
               continue
            if not gene_store.has_key(rank) :
               gene_store[rank]={}
            if not gene_store[rank].has_key(tid) :
               gene_store[rank][tid]=[]
            gene_store[rank][tid].append( (taxid,rc) )


for rank in store.keys() :
   fsname=os.path.basename(fsfile)
   fsfileout=odir+"/"+fsname+"."+rank
   print "create fastsummary file ",fsfileout, "for rank=",rank 
   fh=open(fsfileout,"w")
   save=[]
   for tid in store[rank].keys() :
      if orig.has_key(tid) :
         name_str=orig[tid]
      else :
         name_str=names[tid]
      ### taxid sum
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
      
      gene_lst = []
      rrna_csum=0
      if gene_store.has_key(rank) and gene_store[rank].has_key(tid) :
         gene_lst=gene_store[rank][tid]
         for taxid,count in gene_lst :
            rrna_csum += int(count)
      tup=(wrc_sum,count_sum,tid,name_str,rrna_csum,strain_info)
      save.append(tup)
   sval=sorted(save, key=lambda val : val[0],reverse=True)
   for val in sval :
      if gsfile != "" :
         pcnt = float(val[4])/float(val[1])
         fstr ="%.4f" % pcnt
         out_str=str(val[0])+"\t"+str(val[1]) + "\t" + fstr + "\t" + str(val[2])+"\t"+val[3] + val[5]
      else :
         out_str=str(val[0])+"\t"+str(val[1]) + "\t" + str(val[2])+"\t"+val[3] + val[5]
      fh.write(out_str +"\n")
