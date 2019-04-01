#!/bin/python
import mdtraj as md
import numpy as np

def csd(traj,contact):
# traj=md.load_frame(trajectory,top=structure,index=frame)

 radii={}

 file=open(contact,"r")		### read contact parameters from a file
 for line in file:
  w=line.split()
  if(w[0]=="Ac"):
   Ac=float(w[1])
  elif(w[0]=="Dc"):
   Dc=float(w[1])
  else:
   radii[w[0]]=float(w[1])
 file.close()

 pn=-1		### pair number
 c=[]
 cxyz=[]

 bx=traj.unitcell_lengths[0][0]
 by=traj.unitcell_lengths[0][1]
 bz=traj.unitcell_lengths[0][2]

 for r1 in range(0,traj.n_residues-1):	### molecule A

  x1=traj.xyz[0,r1,0]
  y1=traj.xyz[0,r1,1]
  z1=traj.xyz[0,r1,2]

  for r2 in range(r1+1,traj.n_residues):	### molecule B

   x2=traj.xyz[0,r2,0]
   y2=traj.xyz[0,r2,1]
   z2=traj.xyz[0,r2,2]

   x2=x2-bx*int((x2-x1)/bx)		### deals with unwraped trajectories
   y2=y2-by*int((y2-y1)/by)
   z2=z2-bz*int((z2-z1)/bz)

   if((x2-x1)>(0.5*bx)):		### find the closest image of two
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

#   dxyz=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)	### closest distance
   dxyz2=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2

   rn1=traj.topology.atom(r1).name
   rn2=traj.topology.atom(r2).name

#   rc=np.sqrt((Ac*(radii[rn1]+radii[rn2])*0.5+Dc)**2)		### contact distance criterion
   rc2=(Ac*(radii[rn1]+radii[rn2])*0.5+Dc)**2

#   print(r1,rn1,"-",r2,rn2,rc,dxyz)

   if(dxyz2<rc2):
#    print("fr:",0,"pair:",r1,"-",r2)

    pn=pn+1

    c1=-1
    c2=-1

    if(pn==0):					### first pair
     c.append([r1,r2])
     cxyz.append([[x1,y1,z1],[x2,y2,z2]])
    else:
     for i in range(0,len(c)):			### compare with existing clusters
      if(r1 in c[i]):
       c1=i
      if(r2 in c[i]):
       c2=i
     if((c1==-1)&(c2==-1)):			### doesn't appear in any - add new
      c.append([r1,r2])
      cxyz.append([[x1,y1,z1],[x2,y2,z2]])
     elif((c1!=-1)&(c2!=-1)&(c1!=c2)):		### appears in two different
      c.append(c[c1]+c[c2])
      cxyz.append(cxyz[c1]+cxyz[c2])
      if(c1>c2):
       del c[c1],cxyz[c1]
       del c[c2],cxyz[c2]
      else:
       del c[c2],cxyz[c2]
       del c[c1],cxyz[c1]
     elif((c1!=-1)&(c2==-1)):				### appears in one
      c[c1].append(r2)
      cxyz[c1].append([x2,y2,z2])
     elif((c2!=-1)&(c1==-1)):				### appears in one
      c[c2].append(r1)
      cxyz[c2].append([x1,y1,z1])
     
#    print("clusters:",c)
#    print("cluster sizes:",len(c))

 for m in range(0,traj.n_atoms):			### add monomers
  if(not any(m in i for i in c)):
   c.append([m])
   cxyz.append([traj.xyz[0,m,0],traj.xyz[0,m,1],traj.xyz[0,m,2]])
# print(c)
 return c,cxyz
