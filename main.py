#!/bin/python

import argparse
import cluster
import traj
import genpdb
import csd
import cog
import rdf

parser=argparse.ArgumentParser(description="Parse options.")

parser.add_argument("-s",type=str,help="structure file")
parser.add_argument("-t",type=str,help="trajecotry file")
parser.add_argument("-f",type=int,help="frame number")
parser.add_argument("-c",type=str,help="contact criterion parameters")

args=parser.parse_args()



### reads the structure and a frame from the trajectory
trajectory=traj.traj(args.s,args.t,args.f)

### gives back clustered proteins numbers and their wrapped coordinates
clust,clust_xyz=cluster.cluster(trajectory,args.c)


### determines cluster size distribution function
csd.csd(clust)

### generates a pdb of proteins with wrapped coordinates
genpdb.genpdb(clust,clust_xyz,trajectory)

### calculates centers of geometry of clusters
cog=cog.cog(clust,clust_xyz,trajectory)

### calculates radial distribution function of clusters
rdf=rdf.rdf(clust,clust_xyz,trajectory,cog)

print(clust)
#print(sorted(clust,key=len,reverse=True))
#print(clust.sort(key=len))
