#!/bin/python

import numpy
import argparse
import cluster
import traj
import genpdb
import csd
import cog
import rdf
import os

parser=argparse.ArgumentParser(description="Parse options.")

parser.add_argument("-s",type=str,help="structure file")
parser.add_argument("-t",type=str,help="trajecotry file")
parser.add_argument("-f",type=int,help="frame number")
parser.add_argument("-fndx",type=str,help="file with frame numbers")
parser.add_argument("-c",type=str,help="contact criterion parameters")

args=parser.parse_args()

### determine frame/frames
if(args.fndx!=None):			### file with frame numbers
 d=numpy.loadtxt(str(args.fndx))
elif(args.f!=None):			### frame number
 d=args.f
else:					### default frame number
 d=1

for n in d:
 frame=int(n-1)
 dirout="f"+str(int(n))
 if not os.path.exists(str(dirout)):
  os.mkdir(dirout)

### reads the structure and a frame from the trajectory
 trajectory=traj.traj(args.s,args.t,frame)

### gives back clustered proteins numbers and their wrapped coordinates
 clust,clust_xyz=cluster.cluster(trajectory,args.c)
 cluster.cluster_write(clust,dirout)

### determines cluster size distribution function
 clustdist=csd.csd(clust)
 csd.csd_write(clustdist,dirout)

### generates a pdb of proteins with wrapped coordinates
 genpdb.genpdb(clust,clust_xyz,trajectory,dirout)

### calculates centers of geometry of clusters
 center=cog.cog(clust,clust_xyz,trajectory)
 cog.cog_write(center,dirout)

### calculates radial distribution function of clusters
 rdf.rdf(clust,clust_xyz,trajectory,center,dirout)

#print(clust)
#print(sorted(clust,key=len,reverse=True))
#print(clust.sort(key=len))
