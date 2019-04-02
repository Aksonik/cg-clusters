#!/bin/python

import argparse
import oj
import readtraj
import genpdb
import csdplt
import cog
import rdf

parser=argparse.ArgumentParser(description="Parse options.")

parser.add_argument("-s",type=str,help="structure file")
parser.add_argument("-t",type=str,help="trajecotry file")
parser.add_argument("-f",type=int,help="frame number")
parser.add_argument("-c",type=str,help="contact criterion parameters")

args=parser.parse_args()

trajectory=readtraj.traj(args.s,args.t,args.f)
clust,clust_xyz=oj.csd(trajectory,args.c)

#print(clust,clust_xyz)

csdplt.csdplt(clust)

genpdb.genpdb(clust,clust_xyz,trajectory)

cog=cog.cog(clust,clust_xyz,trajectory)

rdf=rdf.rdf(clust,clust_xyz,trajectory,cog)

