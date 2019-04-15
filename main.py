#!/bin/python

import numpy
import argparse
import cluster
import traj
import genpdb
import csd
import csd_plt
import cog
import rdf
import os

parser=argparse.ArgumentParser(description="Parse options.")

parser.add_argument("-s",type=str,help="structure file (PDB)")
parser.add_argument("-t",type=str,help="trajecotry file (DCD)")
parser.add_argument("-f",type=int,help="frame number (singe integer)")
parser.add_argument("-fn",type=str,help="file with frame numbers (integers)")
parser.add_argument("-c",type=str,help="file with contact criterion parameters")
parser.add_argument("-bs",type=float,help="bin size for RDF analysis")

args=parser.parse_args()

### determine frame/frames
if(args.fn!=None):			### file with frame numbers
 d=numpy.loadtxt(str(args.fn),dtype="int")
elif(args.f!=None):			### frame number
 d=[args.f]
else:					### default frame number
 d=[1]

### bin size for RDF
if(args.bs!=None):
 bs=args.bs
else:
 bs=0.5			### [nm]

for n in d:
 print("frame:",n)
 frame=n-1
 dirout="f"+str(n)
 if not os.path.exists(str(dirout)):
  os.mkdir(dirout)

### reads the structure and a frame from the trajectory
 trajectory=traj.traj(args.s,args.t,frame)

### gives back clustered proteins numbers and their wrapped coordinates
 clust,clust_xyz,contacts=cluster.cluster(trajectory,args.c)
 cluster.cluster_write(clust,dirout)
 cluster.contacts_write(contacts,dirout,args.c)

### determines cluster size distribution function
 clustdist=csd.csd(clust)
 csd.csd_write(clustdist,dirout)

### generates a pdb of proteins with wrapped coordinates
 genpdb.genpdb(clust,clust_xyz,trajectory,dirout)

### calculates centers of geometry of clusters
 center=cog.cog(clust,clust_xyz,trajectory)
 cog.cog_write(center,dirout)

### calculates radial distribution function of clusters
 rdf.rdf(clust,clust_xyz,trajectory,center,dirout,bs)

#print(clust)
#print(sorted(clust,key=len,reverse=True))
#print(clust.sort(key=len))

### average and plot cluster size distributions

csd.csd_avg(d)
#csd_plt.csd_plt()
