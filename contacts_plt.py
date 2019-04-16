#!/bin/python

import matplotlib as mpl
mpl.use('Agg')

import numpy
import matplotlib
from pylab import *

fs=24		### font size (for labels is 3 times larger)
lc="blue"	### font color
lw=1		### line weight (when neighboring labels is bolt)

fi=open("contacts/contacts.dat","r")	### data file (label x, label y, value)

fig=figure(figsize=(12,3))
fig.subplots_adjust(left=0.05,bottom=0.05,right=0.95,top=0.95,hspace=0.05,wspace=0.05)

ax=plt.subplot(1,1,1)

plt.axis('off')

val={}
err={}
labx=[]
laby=[]

for line in fi:
 w=line.split()

 val[w[0],w[1]]=w[2]
 err[w[0],w[1]]=w[4]

 if w[0] not in labx:
  labx.append(w[0])
 if w[1] not in laby:
  laby.append(w[1])

valmat=numpy.zeros((len(labx),len(laby)))
errmat=numpy.zeros((len(labx),len(laby)))

for matx in range(0,len(labx)):
 for maty in range(0,len(laby)):

  if (labx[matx],laby[maty]) in val:
   valmat[matx,maty]=val[labx[matx],laby[maty]]
   errmat[matx,maty]=err[labx[matx],laby[maty]]

### plot label x

for matx in range(0,len(labx)):
 ax.text(matx+2,len(laby)+1,labx[matx],ha="center",va="center",fontsize=fs,weight="bold")

### plot label y

for maty in range(0,len(laby)):
 ax.text(1,len(laby)-maty,laby[maty],ha="center",va="center",fontsize=fs,weight="bold")

### plot value

for matx in range(0,len(labx)):

 if(matx==0):
  linew=3*lw
 else:
  linew=lw

 plot([matx+1.5,matx+1.5],[0,10],color=lc,lw=linew)	### vertical line

 for maty in range(0,len(laby)):

  if(maty==len(laby)-1):
   linew=3*lw
  else:
   linew=lw

  plot([0,10],[maty+1.5,maty+1.5],color=lc,lw=linew)	### horizontal line

  if (labx[matx],laby[maty]) in val:	### leave empty if the matrix is symmetric
   ax.text(matx+2,len(laby)-maty,
   str(valmat[matx,maty])+r" $\pm$ "+str(errmat[matx,maty]),
   ha="center",va="center",fontsize=fs)

xlim(0.5,len(labx)+1.5)
ylim(0.5,len(laby)+1.5)

savefig("contacts/contacts.png")
