#!/bin/python
import numpy
import os
import glob

def rdf(clust,clust_xyz,traj,cog,dirout,bs):

 dd=[]
 rr=[]

 for c in range(0,len(cog)):			### cluster

  cx=cog[c][0]					### center of geometry
  cy=cog[c][1]
  cz=cog[c][2]

  for p in range(0,len(clust_xyz[c])):		### protein

   px=clust_xyz[c][p][0]			### coordinates
   py=clust_xyz[c][p][1]
   pz=clust_xyz[c][p][2]

   r=traj.topology.atom(clust[c][p]).name		### residue
   d=numpy.sqrt((cx-px)**2+(cy-py)**2+(cz-pz)**2)		### distance from the center

   if(p==0):
    dd.append([d])
    rr.append([r])
   else:
    dd[c].append(d)
    rr[c].append(r)

### print(dd)	### distances
### print(rr)

 for c in range(0,len(dd)):		### cluster
  for r in set(rr[c]):			### residue type

   drr=[]

   for d in range(0,len(dd[c])):	### protein

    if(rr[c][d]==r):
     drr.append(dd[c][d])

###     print(c,r,d,dd[c][d])
###   print(drr)			### distances for a given residue type

   bn=[]
   for b in range(0,int(max(dd[c])/bs)+1):
    bn.append(0)

   for dr in range(0,len(drr)):
    dndx=int(drr[dr]/bs)			### corresponding bin index
    bn[dndx]=bn[dndx]+1

###   print(bn)

   f=open(str(dirout)+"/rdf/rdf_"+str(c)+"_"+str(len(clust[c]))+"_"+str(len(drr))+"_"+str(r)+".dat","w")

   for b in range(0,len(bn)):
    vol=4.0/3.0*3.141592653589793*((b*bs+bs)**3-(b*bs)**3)	### [nm^3]
    print(b*bs+0.5*bs,bn[b]/vol,file=f)				### [1/nm^3]

   f.close()



def rdf_avg(d,molecules_types):

 for mt in molecules_types:

  ax=[]
  ay=[]

  for i in range(0,1000):
   ax.append([])
   ay.append([])

  for f in d:
   finame=glob.glob("f"+str(f)+"/rdf/rdf_0_5??_*_"+mt+".dat")
   if(os.path.isfile(finame)):

    fi=open(finame[0],"r")

    data=numpy.loadtxt(fi)
    x=data[:,0]
    y=data[:,1]

    for i in range(0,len(x)):
     ax[i].append(x[i])
     ay[i].append(y[i])

    fi.close()

  dd="rdf"
  if not os.path.exists(dd):
   os.mkdir(dd)

  fi=open(dd+"/rdf_"+mt+".dat","w")

  for i in range(0,len(ax)):
   if(ax[i]!=[]):
    print("%9.3f %9.3f %9.3f %9.3f" % (numpy.mean(ax[i]),numpy.mean(ay[i]),
numpy.std(ay[i]),numpy.std(ay[i])/numpy.sqrt(len(ax[i]))),file=fi)

  fi.close()
