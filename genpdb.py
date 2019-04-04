#!/bin/python

def genpdb(clust,cxyz,traj):

 f=open("genpdb.pdb","w")

 bx=traj.unitcell_lengths[0][0]*10.0
 by=traj.unitcell_lengths[0][1]*10.0
 bz=traj.unitcell_lengths[0][2]*10.0

 print("%s%9.3f%9.3f%9.3f%s" % 
 ("CRYST1",bx,by,bz,"  90.00  90.00  90.00 P 1           1"),file=f)

 for c in clust:

  for p in c:
   n=p+1
   x=cxyz[clust.index(c)][c.index(p)][0]*10.0	### [nm] -> [A]
   y=cxyz[clust.index(c)][c.index(p)][1]*10.0
   z=cxyz[clust.index(c)][c.index(p)][2]*10.0
   r=traj.topology.atom(p).name

   print("%6s%5i%5s%4s%6i    %8.3f%8.3f%8.3f%s%s" % 
   ("ATOM  ",n,r,r,n,x,y,z,"  1.00  1.00      ","PRO0"),file=f)

 f.close()
