#!/bin/python

import argparse
import oj
import readtraj
#import rdf

parser=argparse.ArgumentParser(description="Parse options.")

parser.add_argument("-s",type=str,help="structure file")
parser.add_argument("-t",type=str,help="trajecotry file")
parser.add_argument("-f",type=int,help="frame number")
parser.add_argument("-c",type=str,help="contact criterion parameters")

args=parser.parse_args()

trajectory=readtraj.traj(args.s,args.t,args.f)
clusters=oj.csd(trajectory,args.c)
#rdf=rdf.rdf(trajectory,clusters)

print(clusters)
