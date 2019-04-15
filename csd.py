#!/bin/python
import numpy
import os

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



def csd_avg(d):

 ax=[]
 ay=[]

 for i in range(0,10000):
  ax.append([])
  ay.append([])

 for f in d:
  f=open("f"+str(f)+"/csd.dat","r")

  data=numpy.loadtxt(f)
  x=data[:,0]
  y=data[:,1]

  for i in range(0,len(x)):
   ax[i].append(x[i])
   ay[i].append(y[i])

  f.close()

 d="csd"
 if not os.path.exists(d):
  os.mkdir(d)

 f=open(d+"/csd.dat","w")

 for i in range(0,len(ax)):
  if(ax[i]!=[]):
   print("%6i %9.3f %9.3f %9.3f" % (numpy.mean(ax[i]),numpy.mean(ay[i]),
numpy.std(ay[i]),numpy.std(ay[i])/numpy.sqrt(len(ax[i]))),file=f)

 f.close()
