#!/usr/bin/env python

import sys
from operator import itemgetter

summfile=sys.argv[1]
rankfile=sys.argv[2]
#summtype=sys.argv[3] ## strain,species,species
#fastsumm=sys.argv[2]

def loadRank(rfile) :
   save={}
   fh1=open(rfile) 
   for ln in fh1 :
      ln=ln.rstrip()
      vl=ln.split()
      save.setdefault(vl[0],vl[1])
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
   return arr[iter],arr[iter+1],int(arr[iter+4]),int(arr[iter+6]),float(arr[iter+7])

rankMap=loadRank(rankfile)
fh=open(summfile)

names={}
parent={}
child={}
wrdcnt,kmercnt,rdcnt={},{},{}
parent.setdefault(1,1)
lines=[('1',0)]
last_tab_cnt=-1
for line in fh :
   line=line.rstrip()
   vals=line.split('\t')
   if vals[0]=='Name' :
      continue
   num_tabs=cntTabs(vals)
   pn,cnode,val1,val2,val3 = getNode(vals)
   names[cnode] = pn
   kmercnt[cnode]=val1
   rdcnt[cnode]=val2
   wrdcnt[cnode]=val3
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


save_calls=[]
open=['1']
done={}
seen={}
visitor={}
plasmids={}
while open != []:
   cnode=open[0]
   lst=[]
   if child.has_key(cnode) :
      lst=child[ cnode ]
   isDoneChild=True
   print "cnode=",cnode," children=",lst
   it=0
   while it < len(lst) :
      if not seen.has_key(lst[it]) :
         open.insert(0,lst[it])
         seen.setdefault(lst[it],1)
      if not done.has_key(lst[it]) :
         isDoneChild=False
      it += 1
   if isDoneChild :
      open.pop(0)
      done.setdefault(cnode)
      ## process this node
      if int(cnode) >= 10000000 : 
         pval=parent[cnode]
         if rankMap[pval] == "strain" :
            pval = parent[pval]
         if not plasmids.has_key(pval) : 
            plasmids.setdefault(pval,[])
         plasmids[pval].append(cnode)
      elif rankMap[cnode] == "species" :
         tot_wrdcnt = wrdcnt[cnode]
         tot_rdcnt = rdcnt[cnode]
         tot_kmercnt = kmercnt[cnode]
         if not child.has_key(cnode ) : 
            the_call=cnode
         else :
            lst=child[ cnode ]
            slst = sorted( lst, key=lambda val : wrdcnt[val], reverse=True )
            for it in range(0,len(slst),1) : 
               alt=slst[it]
               tot_wrdcnt += wrdcnt[alt]
               tot_rdcnt += rdcnt[alt]
               tot_kmercnt += kmercnt[alt]
               the_call=slst[0]

         print "huh:",lst
         print "hello",the_call,tot_wrdcnt,tot_rdcnt,tot_kmercnt
         save_calls.append( (the_call,cnode,tot_wrdcnt,tot_rdcnt,tot_kmercnt) )
      else :
         print "ignore?",cnode,rankMap[cnode]

rep=sorted ( save_calls, key = itemgetter(2), reverse=True )
for val in rep :
   spec=val[1]
   print names[val[0]],val
   if plasmids.has_key(spec) :
      for val in plasmids[spec] :
         print "plamid:",names[val[0]],val,wrdcnt[spec],rdcnt[spec],kmercnt[spec]
