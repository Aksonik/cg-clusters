#!/bin/python

def csdplt(clust):

 maxc=1			### the largest cluster found
 totp=0			### total number of proteins

 for c in clust:
  totp=totp+len(c)
  if(len(c)>maxc):
   maxc=len(c)

 f=open("csd.dat","w")

 for s in range(1,maxc+1):
  cp=0
  for c in clust:
   if(len(c)==s):
    cp=cp+len(c)
  print(s,cp/totp,file=f)

 f.close()
