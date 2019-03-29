#!/bin/python
import mdtraj as md
import numpy as np

traj=md.load('../cg/sys.dcd',top='../cg/sys.pdb')

pn=-1	### pair number
c=[]

frame=0

bx=traj.unitcell_lengths[frame][0]
by=traj.unitcell_lengths[frame][1]
bz=traj.unitcell_lengths[frame][2]

for r1 in range(0,traj.n_residues-1):	### loops over all protein pairs
 for r2 in range(r1+1,traj.n_residues):

  x1=traj.xyz[frame,r1,0]
  y1=traj.xyz[frame,r1,1]
  z1=traj.xyz[frame,r1,2]

  x2=traj.xyz[frame,r2,0]
  y2=traj.xyz[frame,r2,1]
  z2=traj.xyz[frame,r2,2]

  dx=(x1-x2)%bx				### deals with unwraped trajectories too
  dy=(y1-y2)%by
  dz=(z1-z2)%bz

  if(dx>0.5*bx):
   dx=dx-bx

  if(dy>0.5*by):
   dy=dy-by

  if(dz>0.5*bz):
   dz=dz-bz

#  dxyz=np.sqrt(dx**2+dy**2+dz**2)	### closest distance
  dxyz2=dx**2+dy**2+dz**2

  rn1=traj.topology.atom(r1).name
  rn2=traj.topology.atom(r2).name

  if(rn1=="CGA"):
   rc1=1.81273
  elif(rn1=="CGB"):
   rc1=2.06692
  elif(rn1=="CGC"):
   rc1=2.33446
  elif(rn1=="CGP"):
   rc1=2.33446
  elif(rn1=="CGL"):
   rc1=2.69569

  if(rn2=="CGA"):
   rc2=1.81273
  elif(rn2=="CGB"):
   rc2=2.06692
  elif(rn2=="CGC"):
   rc2=2.33446
  elif(rn2=="CGP"):
   rc2=2.33446
  elif(rn2=="CGL"):
   rc2=2.69569

  Ac=1.0
  Dc=0.7

#  rc=np.sqrt((Ac*(rc1+rc2)*0.5+Dc)**2)		### contact criterion
  rc2=(Ac*(rc1+rc2)*0.5+Dc)**2

#  print(r1+1,rn1,"-",r2+1,rn2,rc,dxyz)

  if(dxyz2<rc2):
   print("fr:",frame,"pair:",r1,"-",r2)

   pn=pn+1

   c1=-1
   c2=-1

   if(pn==0):					### first pair
    c.append([r1,r2])
   else:
    for i in range(0,len(c)):			### compare with existing clusters
     if(r1 in c[i]):
      c1=i
     if(r2 in c[i]):
      c2=i
    if((c1==-1)&(c2==-1)):			### doesn't appear in any - add new
     c.append([r1,r2])
    elif((c1!=-1)&(c2!=-1)&(c1!=c2)):		### appears in two different
     c.append(c[c1]+c[c2])
     if(c1>c2):
      del c[c1]
      del c[c2]
     else:
      del c[c2]
      del c[c1]
    elif((c1!=-1)&(c2==-1)):				### appears in one
     c[c1].append(r2)
    elif((c2!=-1)&(c1==-1)):				### appears in one
     c[c2].append(r1)
     
   print("clusters:",c)
   print("cluster sizes:",len(c))

for m in range(0,traj.n_atoms):			### add monomers
 if(not any(m in i for i in c)):
  c.append([m])

print(c)

