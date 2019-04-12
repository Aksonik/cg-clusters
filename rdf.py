#!/bin/python
import numpy as np
import os

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
   d=np.sqrt((cx-px)**2+(cy-py)**2+(cz-pz)**2)		### distance from the center

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
