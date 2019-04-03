#!/bin/python
import mdtraj as md
import numpy as np

def cluster(traj,contact):

 file=open(contact,"r")		### read contact parameters from a file

 radii={}

 for line in file:
  w=line.split()

  if(w[0]=="Ac"):
   Ac=float(w[1])
  elif(w[0]=="Dc"):
   Dc=float(w[1])
  else:
   radii[w[0]]=float(w[1])

 file.close()

 c=[]		### clustered proteins
 cxyz=[]	### coordinates of clustered proteins

 bx=traj.unitcell_lengths[0][0]
 by=traj.unitcell_lengths[0][1]
 bz=traj.unitcell_lengths[0][2]

 pr=[]					### remaining proteins
 for p in range(0,traj.n_residues):
  pr.append(p)

 cn=-1

 for p0 in pr:	### initial protein (one for a cluster)

  cn=cn+1

  x0=traj.xyz[0,p0,0]
  y0=traj.xyz[0,p0,1]
  z0=traj.xyz[0,p0,2]

  c.append([p0])
  cxyz.append([[x0,y0,z0]])

  pr.remove(p0)

  for p1 in c[cn]:	### proteins in a cluster

   p1ndx=c[cn].index(p1)

   x1=cxyz[cn][p1ndx][0]
   y1=cxyz[cn][p1ndx][1]
   z1=cxyz[cn][p1ndx][2]

   rn1=traj.topology.atom(p1).name

   for p2 in pr:			### proteins not in clusters

    x2=traj.xyz[0,p2,0]
    y2=traj.xyz[0,p2,1]
    z2=traj.xyz[0,p2,2]

    rn2=traj.topology.atom(p2).name

    x2=x2-bx*int((x2-x1)/bx)		### deals with unwraped trajectories
    y2=y2-by*int((y2-y1)/by)
    z2=z2-bz*int((z2-z1)/bz)

    if((x2-x1)>(0.5*bx)):		### find the closest image of the two
     x2=x2-bx
    elif((x2-x1)<(-0.5*bx)):
     x2=x2+bx

    if((y2-y1)>(0.5*by)):
     y2=y2-by
    elif((y2-y1)<(-0.5*by)):
     y2=y2+by

    if((z2-z1)>(0.5*bz)):
     z2=z2-bz
    elif((z2-z1)<(-0.5*bz)):
     z2=z2+bz

#    dxyz=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)	### closest distance
    dxyz2=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2

#    rc=np.sqrt((Ac*(radii[rn1]+radii[rn2])*0.5+Dc)**2)		### contact distance criterion
    rc2=(Ac*(radii[rn1]+radii[rn2])*0.5+Dc)**2

#    print("%2i-%2i (%3s-%3s) %12.6f %12.6f" % (p1,p2,rn1,rn2,rc,dxyz))
#    print(c,pr)

    if(dxyz2<rc2):
     c[cn].append(p2)
     cxyz[cn].append([x2,y2,z2])
     pr.remove(p2)

     print("proteins left:",len(pr))

 for m in pr:			### add monomers
  c.append([m])
  cxyz.append([[traj.xyz[0,m,0],traj.xyz[0,m,1],traj.xyz[0,m,2]]])
  pr.remove(m)

  print("proteins left:",len(pr))
 
 return(c,cxyz)
