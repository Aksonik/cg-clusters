#!/bin/python

import matplotlib as mpl
mpl.use('Agg')

import numpy as np
from pylab import *

def rdf_plot(molecules_types):

 colors=["blue","red","green"]

 fig=figure(figsize=(10,5))
 fig.subplots_adjust(left=0.13,bottom=0.16,right=0.97,top=0.95,hspace=0.05,wspace=0.05)

 c=-1

 for mt in molecules_types:

  c=c+1

  f="rdf/rdf_"+str(mt)+".dat"

  d=np.loadtxt(f)
  x=d[:,0]
  y=d[:,1]
  yerr=d[:,3]

  xx=np.insert(x,0,x[0]-1)		### to have first step full
  yy=np.insert(y,0,y[0])

  b=x[1]-x[0]

  plt.step(xx+b*0.5,yy,color=colors[c],linewidth=2.0,label=str(mt))
  plt.errorbar(x,y,yerr=yerr,ecolor=colors[c],ls="none",elinewidth=1.0,capsize=4)
  plt.plot(x,y,marker="s",color=colors[c],ls="None")

 plt.xlabel(r'radius [nm]',fontsize=24)
 plt.ylabel(r'density [1/nm$^3$]', fontsize=24)

# plt.xlim(-1,25)
# plt.ylim(0,20)
# plt.xticks(np.arange(0,550,100),fontsize=20)
# plt.yticks(np.arange(0,20,2),fontsize=20)

 legend(loc=1,fontsize=20,fancybox=True).get_frame().set_alpha(0.5)

 plt.xticks(fontsize=20)
 plt.yticks(fontsize=20)

 plt.savefig("rdf/rdf.png")
