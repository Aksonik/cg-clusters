#!/bin/python

import numpy
import argparse
import cluster
import traj
import genpdb
import csd
import cog
import rdf
import rdf_plt
import os
import solublim
import cluster_traceback

### parameters

parser=argparse.ArgumentParser(description="Parse options.")

parser.add_argument("-s",type=str,help="structure file (PDB)")
parser.add_argument("-t",type=str,help="trajecotry file (DCD)")
parser.add_argument("-f",type=int,help="frame number (singe integer)")
parser.add_argument("-fn",type=str,help="file with frame numbers (integers)")
parser.add_argument("-c",type=str,help="file with contact criterion parameters")
parser.add_argument("-bs",type=float,help="bin size for RDF analysis")
parser.add_argument("-sl",type=int,help="cluster size threshold for the solubility limit")

args=parser.parse_args()

### read frame/frames for the analysis

if(args.fn!=None):			### file with frame numbers
 d=numpy.loadtxt(str(args.fn),dtype="int")
elif(args.f!=None):			### frame number
 d=[args.f]
else:
 d=[1]					### default frame number

### bin size for radial distribution function

if(args.bs!=None):
 bs=args.bs
else:
 bs=0.5					### default size [nm]

### loop over the frames

for n in d:
 print("frame:",n)
 frame=n-1
 dirout="f"+str(n)
 if not os.path.exists(str(dirout)):
  os.mkdir(dirout)

### read a structure (PDB) and a frame (coordinates) from the trajectory (DCD)

 trajectory=traj.traj(args.s,args.t,frame,args.c)

### read the box size [nm]

 bx,by,bz=cluster.box(trajectory,dirout)

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

### calculates the solubility limit

 if(args.sl!=None):
  solublim.solublim(dirout,args.c,args.t,args.s,args.sl)





#print(clust)
#print(sorted(clust,key=len,reverse=True))
#print(clust.sort(key=len))

### average and plot cluster size distributions

csd.csd_avg(d)

import csd_plt
csd_plt

### average and plot number of contacts

cluster.contacts_avg(d)
import contacts_plt
contacts_plt

### average and plot radial distrinution function
### of which cluster? - needs to be done

#mt=cluster.molecules_types(args.c)
#rdf.rdf_avg(d,mt)
#rdf_plt.rdf_plot(mt)


### loop over the frames

for n in d:
 print("frame:",n)
 frame=n-1
 dirout="f"+str(n)
 if not os.path.exists(str(dirout)):
  os.mkdir(dirout)

### trace back a cluster

 cluster_traceback.cluster_traceback()
 cluster_traceback.cluster_traceback_write(n)


