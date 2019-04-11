#!/bin/python

def csd(clust):

 maxc=1			### the largest cluster size found
 totp=0			### total number of proteins

 for c in clust:
  totp=totp+len(c)
  if(len(c)>maxc):
   maxc=len(c)

 csd=[]

 for s in range(1,maxc+1):	### luster size
  cp=0
  for c in clust:
   if(len(c)==s):
    cp=cp+len(c)

  csd.append([s,cp/totp*100])	### [%]

 return csd



def csd_write(csd,dirout):

 f=open(str(dirout)+"/csd.dat","w")

 for c in csd:
  print(c[0],c[1],file=f)

 f.close()
