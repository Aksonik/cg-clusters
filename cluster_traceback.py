#!/bin/python

import sys



### a cluster that is going to be traced back

ff=6000
ll=0

fi=open("f"+str(ff)+"/cluster.dat","r")

l=-1
for line in fi:
 l=l+1
 if(l==0):
  w=line.split(" ")
  c=w[1:len(w)-1]	### not the first (size) and not the last (new line)
  cc=c
  break

fi.close()



### compare with clusters from the previous frames

for f in range(6000,0,-50):
 fi=open("f"+str(f)+"/cluster.dat","r")
 fo=open("f"+str(f)+"/cluster_traceback.dat","w")

 permax=0

 for line in fi:

  w=line.split(" ")
  c=w[1:len(w)-1]	### not the first (size) and not the last (new line)

  co=(set(c)&set(cc))	### molecules common for the previous and current clusters
  per=float(len(co))/float(len(cc))

  if(per>permax):
   permax=per
   cmax=c
   ccmax=c

# print("previous cluster:",cc)
# print("current cluster:",ccmax)
# print("common molecules: %12.6f" % permax)

 cc=ccmax

 print("%9i %6.3f %s" % (f,permax,cc),file=fo)

 fi.close()
 fo.close()
