#!/bin/python
import mdtraj as md
import cluster

def traj(structure,trajectory,frame,contact):

### in the trajecotry take into account only atoms 
### that are defined in the contact criterion file
### i.e. contact.dat

 molecule_types=cluster.molecules_types(contact)

 fc=open(structure,"r")

 atomsndx=[]

 i=-1
 for line in fc:
  w=line.split()

  if(w[0]=="ATOM"):
   i=i+1
   if w[2] in molecule_types:
    atomsndx.append(i)

 fc.close()

### topsel=md.load(str(structure)).topology
### atomsndx=topsel.select('resname CGA or resname CGB or resname CGC')

 traj=md.load_frame(trajectory,top=structure,index=frame,atom_indices=atomsndx)
 return traj
