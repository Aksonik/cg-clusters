#!/bin/python
import mdtraj as md
import numpy as np

### xyz of two atoms
def dist(c1,c2):
 return np.sqrt(np.sum((c1-c2)**2))

traj=md.load('../cg/sys.dcd',top='../cg/sys.pdb')
#traj=md.load('../aa/sys.dcd',top='../aa/sys.pdb')

print(traj)

print('How many atoms?    %s' % traj.n_atoms)
print('How many residues? %s' % traj.n_residues)
print('What is the box size in x? %s' % traj.unitcell_lengths[0][0])	### [nm]

frame_ndx=2 	### from zero
atom_ndx=10	### from zero

print('x: %s\ty: %s\tz: %s' % tuple(traj.xyz[frame_ndx, atom_ndx,:]))	### [nm]

topology = traj.topology
print(topology)



print('Fifth atom: %s' % topology.atom(4))
print('All atoms: %s' % [atom for atom in topology.atoms])

print('Second residue: %s' % traj.topology.residue(1))
print('All residues: %s' % [residue for residue in traj.topology.residues])

atom = topology.atom(10)
print('''Hi! I am the %sth atom, and my name is %s. 
I am a %s atom with %s bonds. 
I am part of an %s residue.''' % ( atom.index, atom.name, atom.element.name, atom.n_bonds, atom.residue.name))

print([atom.index for atom in topology.atoms if atom.element.symbol is 'CGB' and atom.is_sidechain])

print([residue for residue in topology.chain(0).residues if residue.index % 2 == 0])


print(topology.select('resid 1 to 2'))
print(topology.select('name N and backbone'))

selection = topology.select_expression('name CGB and resid 1 to 5')
print(selection)

###############################

cn=-1
c=[]

frame=0

bx=traj.unitcell_lengths[frame][0]
by=traj.unitcell_lengths[frame][1]
bz=traj.unitcell_lengths[frame][2]

for r1 in range(0,traj.n_residues-1):	### loops over all protein pairs
 for r2 in range(r1+1,traj.n_residues):

  x1=traj.xyz[frame,r1,0]
  y1=traj.xyz[frame,r1,1]
  z1=traj.xyz[frame,r1,2]

  x2=traj.xyz[frame,r2,0]
  y2=traj.xyz[frame,r2,1]
  z2=traj.xyz[frame,r2,2]

  dx=(x1-x2)%bx				### deals with unwraped trajectories too
  dy=(y1-y2)%by
  dz=(z1-z2)%bz

  if(dx>0.5*bx):
   dx=dx-bx

  if(dy>0.5*by):
   dy=dy-by

  if(dz>0.5*bz):
   dz=dz-bz

#  dxyz=np.sqrt(dx**2+dy**2+dz**2)	### closest distance
  dxyz2=dx**2+dy**2+dz**2

  rn1=topology.atom(r1).name
  rn2=topology.atom(r2).name

  if(rn1=="CGA"):
   rc1=1.81273
  elif(rn1=="CGB"):
   rc1=2.06692
  elif(rn1=="CGC"):
   rc1=2.33446
  elif(rn1=="CGP"):
   rc1=2.33446
  elif(rn1=="CGL"):
   rc1=2.69569

  if(rn2=="CGA"):
   rc2=1.81273
  elif(rn2=="CGB"):
   rc2=2.06692
  elif(rn2=="CGC"):
   rc2=2.33446
  elif(rn2=="CGP"):
   rc2=2.33446
  elif(rn2=="CGL"):
   rc2=2.69569

  Ac=1.0
  Dc=0.7

#  rc=np.sqrt((Ac*(rc1+rc2)*0.5+Dc)**2)		### contact criterion
  rc2=(Ac*(rc1+rc2)*0.5+Dc)**2

#  print(r1+1,rn1,"-",r2+1,rn2,rc,dxyz)

  if(dxyz2<rc2):
   print(frame+1,r1+1,r2+1)

   cn=cn+1
   c.append([cn,r1])
   c.append([cn,r2])

print(c)

c=[]

#print("distance:",np.sqrt(np.sum((traj.xyz[frame,r1,:]-traj.xyz[frame,r2,:])**2)))
#print("distance:",dist(traj.xyz[frame,r1,:],traj.xyz[frame,r2,:]))
