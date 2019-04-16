#!/bin/python

import matplotlib as mpl
mpl.use('Agg')

import numpy as np
from pylab import *

colors=["blue"]

fig=figure(figsize=(16,4))
fig.subplots_adjust(left=0.06,bottom=0.18,right=0.97,top=0.95,hspace=0.05,wspace=0.05)

f="csd/csd.dat"
d=np.loadtxt(f)
x=d[:,0]
y=d[:,1]
yerr=d[:,3]

xx=np.insert(x,0,x[0]-1)		### to have first step full
yy=np.insert(y,0,y[0])

plt.step(xx+0.5,yy,color=colors[0],linewidth=2.0)
plt.errorbar(x,y,yerr=yerr,ecolor=colors[0],ls="none",elinewidth=1.0)
plt.plot(x,y,marker="s",color=colors[0],ls="None")

plt.xlabel(r'cluster size',fontsize=24)
plt.ylabel(r'proteins [%]', fontsize=24)

plt.xlim(-10,550)
plt.ylim(0,20)
plt.xticks(np.arange(0,550,100),fontsize=20)
plt.yticks(np.arange(0,20,2),fontsize=20)

#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)

plt.savefig("csd/csd.png")
