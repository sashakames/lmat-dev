#!/usr/bin/env python

import sys
from operator import itemgetter

summfile=sys.argv[1]
rankfile=sys.argv[2]
fsummfile=sys.argv[3]
plasmidFile=sys.argv[4]
out_base=sys.argv[5]
rank_calls=sys.argv[6]

globalPlasmids={}
def loadPlasmid(file) :
   fh1=open(file) 
   for ln in fh1 :
      ln=ln.rstrip()
      globalPlasmids.setdefault(int(ln),1)
   #return save

loadPlasmid(plasmidFile)

def isPlasmid(id) :
   if id >= 10000000 or globalPlasmids.has_key(id) :
      return True
   return False

def loadRank(rfile) :
   save={}
   fh1=open(rfile) 
   for ln in fh1 :
      ln=ln.rstrip()
      vl=ln.split()
      save.setdefault(int(vl[0]),vl[1])
   return save

def loadFastSumm(rfile) :
   save={}
   fh1=open(rfile)
   for ln in fh1 :
      ln=ln.rstrip()
      vl=ln.split()
      save.setdefault(int(vl[2]),ln)
   return save

def cntTabs(arr) :
   cnt=0
   for iter in range(len(arr)) :
      if arr[iter] != '' :
         break
      cnt+=1 
   return cnt

def getNode(arr) :
   cnt,iter=0,0
   #print "check",arr

   while iter < len(arr) :
      if arr[iter] != '' :
         break
      iter += 1
   #print "stop",iter,"[",arr[iter],"]"
   return arr[iter],int(arr[iter+1]),int(arr[iter+4]),int(arr[iter+6]),float(arr[iter+7]),int(arr[iter+5])

rankMap=loadRank(rankfile)
fsumMap=loadFastSumm(fsummfile)
fh=open(summfile)

names={}
parent={}
child={}
wrdcnt,kmercnt,rdcnt,call_kmercnt={},{},{},{}
parent.setdefault(1,1)
lines=[('1',0)]
last_tab_cnt=-1
for line in fh :
   line=line.rstrip()
   vals=line.split('\t')
   if vals[0]=='Name' :
      continue
   num_tabs=cntTabs(vals)
   pn,cnode,val1,val2,val3,val4 = getNode(vals)
   names[cnode] = pn
   kmercnt[cnode]=val1
   rdcnt[cnode]=val2
   wrdcnt[cnode]=val3
   call_kmercnt[cnode]=val4
   ## find parent for current node
   while lines != [] : 
      pnode,last_tab_cnt=lines[0]
      #print "tabcnt",pnode,last_tab_cnt, cnode, num_tabs
      if num_tabs > last_tab_cnt :
         if not child.has_key(pnode) :
            child.setdefault(pnode,[])
         child[ pnode ].append( cnode )
         parent[ cnode ] = pnode
         #print "setting ",cnode," parent is", pnode
         #print "children ",pnode," is", child[pnode]
         break
      else :
         #print "should back up the tree"
         lines.pop(0)
   lines.insert(0,(cnode,num_tabs))



def summNode(cnode,call_type) :
   if (rankMap.has_key(cnode) and rankMap[cnode] == call_type and not isPlasmid(cnode)) or (isPlasmid(cnode) and call_type == "plasmid") :
      if cnode == 10000179 :
         print "huh?",cnode
      tot_wrdcnt = wrdcnt[cnode]
      tot_rdcnt = rdcnt[cnode]
      tot_kmercnt = kmercnt[cnode]
      the_call=cnode

      save_strains = []
      lopen=[]
      if child.has_key(cnode) :
         lopen=child[ cnode ]
      while lopen != [] :
         alt = lopen.pop(0)
         ## for species calls don't use plasmid reads - report these separately 
         if (call_type == "species" and not isPlasmid(alt)) or (call_type != "species" ) :
            tot_wrdcnt += wrdcnt[alt]
            tot_rdcnt += rdcnt[alt]
            tot_kmercnt += kmercnt[alt]
         if call_type == "species" and rankMap.has_key(alt) and rankMap[alt] == "strain" and not isPlasmid(cnode):
            save_strains.append( alt )

         if child.has_key(alt) :
            for nd in child[alt] :
               lopen.append(nd)
         
      if save_strains != [] :
         slst = sorted( save_strains, key=lambda val : wrdcnt[val], reverse=True )
         ## for now this returns just the top strain call
         the_call=slst[0]
      return (cnode,the_call,tot_wrdcnt,tot_rdcnt,tot_kmercnt)
      
def bread_first_traverse(call_type) :
   save_calls=[]
   open=[1]
   while open != []:
      cnode=open.pop(0)
      ## some plasmids are ranked as "species", still skip these during species call
      if (call_type == "plasmid" and isPlasmid(cnode)) or (call_type == rankMap[cnode] and not isPlasmid(cnode)) :
         #if cnode >= 10000000 :
            #print "huh",cnode,rankMap[cnode]
            #sys.exit(0)
         result=summNode(cnode,call_type)
         save_calls.append(result)
      else :
         ## should be higher order
         lst=[]
         if child.has_key(cnode) :
            lst=child[ cnode ]
         for nd in lst :
            open.insert(0,nd)
   return save_calls


def doPrn(save_calls,outh) :
   rep=sorted ( save_calls, key = itemgetter(2), reverse=True )
   pstr="% of Reads, Avg Read Score, Weighted Read Count (WRC), Read Count (RC), k-mer count, k-mer %, Original WRC, Original RC, Name, Taxid"
   outh.write(pstr+"\n")
   pstr=str(tot_pcnt)+"\t"+str(avg)+"\t"+str(wrc)+"\t"+str(rc)+"\t"+str(kmrcnt)+"\t"+str(ck)+"\t"+str(kpcnt)+"\t"+str(owrc)+"\t"+str(orc)+"\t"+prnName+"\t"+str(call_id)

   rc_sum=0
   for val in rep :
      rc_sum+=val[3]
   for val in rep :
      call_id=val[1]
      owrc,orc=-1,-1
      if fsumMap.has_key(call_id) : 
         v1=fsumMap[call_id].split('\t')
         prnName=v1[3]
         owrc,orc=v1[0],v1[1]
      else :
         prnName=names[call_id]
      wrc,rc,kmrcnt=val[2],val[3],val[4]
      avg=float(wrc)/float(rc)
      ck=call_kmercnt[call_id]
      kpcnt=float(kmrcnt)/float(ck)
      tot_pcnt=float(rc)/float(rc_sum)
      pstr=str(tot_pcnt)+"\t"+str(avg)+"\t"+str(wrc)+"\t"+str(rc)+"\t"+str(kmrcnt)+"\t"+str(ck)+"\t"+str(kpcnt)+"\t"+str(owrc)+"\t"+str(orc)+"\t"+prnName+"\t"+str(call_id)
      outh.write(pstr+"\n")

for ranktype in rank_calls.split() :
   outfile=open(out_base+"."+ranktype,"w")
   save_calls=bread_first_traverse(ranktype)
   doPrn(save_calls,outfile)
