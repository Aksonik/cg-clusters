#!/bin/python
import numpy
import cluster
import mdtraj

def solublim(dirout,contact,trajectory,structure,cst):

### read the box size and calculate the volume

 fi=open(str(dirout)+"/box.dat","r")
 data=numpy.loadtxt(fi)
 bx=data[0]
 by=data[1]
 bz=data[2]
 fi.close()

 bv=bx*by*bz

### read molecules radii

 radii=cluster.molecules_radii(contact)

### read and divide molecules between large and small clusters

 lc=[]		### molecules in large clusters
 sc=[]		### molecules in small clusters

 fi=open(str(dirout)+"/cluster.dat","r")
 for line in fi:
  w=line.split()
  c=w[1:]		### skip the first item, i.e. cluster size
  if(len(c)>cst):
   for j in c:
    lc.append(j)
  else:
   for j in c:
    sc.append(j)
 fi.close()

### read molecule types and calculate the volume

 lv=0.0

 traj=mdtraj.load_frame(trajectory,top=structure,index=0)
 for lm in lc:
  res=traj.topology.atom(int(lm)).name
  rad=radii[res]

  lv+=4.0/3.0*3.141592653589793*rad**3

### solubility limit [mM]
### number of proteins in small clusters per saturated volume, 
### i.e. box volume minus volume of proteins in large clusters

### 1 prot. / 1 nm3 = 1.660539 [M]

 fi=open(str(dirout)+"/solublim.dat","w")
 print(len(sc)/(bv-lv)*1.660539*1000,file=fi)
 fi.close()
