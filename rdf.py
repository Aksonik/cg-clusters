#!/bin/python

def rdf(traj,clust):

 bx=traj.unitcell_lengths[0][0]
 by=traj.unitcell_lengths[0][1]
 bz=traj.unitcell_lengths[0][2]

 c=[]

 for cs in range(0,len(clust)):		### loop over clusters
  pn=-1
  for cn in clust[cs]:			### loop over proteins
   pn=pn+1

   x=traj.xyz[0,cn,0]
   y=traj.xyz[0,cn,1]
   z=traj.xyz[0,cn,2]
   
   if(pn==0):				### add first protein
    c.append([[x,y,z]])
   else:
    for pc in range(0,len(c)):			### compare added proteins
     for pq in range(pc,len(clust[cs])):	### with the rest

      for nx in [-1,0,1]:
       for ny in [-1,0,1]:
        for nz in [-1,0,1]:
         


#    c[cs].append([x1,y1,z1])
 print(c)

"""

   x2=traj.xyz[0,r2,0]
   y2=traj.xyz[0,r2,1]
   z2=traj.xyz[0,r2,2]

   dx=(x1-x2)%bx				### deals with unwraped trajectories too
   dy=(y1-y2)%by
   dz=(z1-z2)%bz

   if(dx>0.5*bx):
    dx=dx-bx

   if(dy>0.5*by):
    dy=dy-by

   if(dz>0.5*bz):
    dz=dz-bz

#   dxyz=np.sqrt(dx**2+dy**2+dz**2)	### closest distance
   dxyz2=dx**2+dy**2+dz**2

   rn1=traj.topology.atom(r1).name
   rn2=traj.topology.atom(r2).name
"""
