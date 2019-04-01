#!/bin/python
import mdtraj as md

def traj(structure,trajectory,frame):
 traj=md.load_frame(trajectory,top=structure,index=frame)
 return traj
