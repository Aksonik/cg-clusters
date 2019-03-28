#!/bin/python
import mdtraj as md

traj=md.load('sys.dcd',top='sys.pdb')

print(traj)


print('How many atoms?    %s' % traj.n_atoms)
print('How many residues? %s' % traj.n_residues)

frame_ndx=2 	### from zero
atom_ndx=10	### from zero

print('x: %s\ty: %s\tz: %s' % tuple(traj.xyz[frame_ndx, atom_ndx,:]))	### [nm]

topology = traj.topology
print(topology)



print('Fifth atom: %s' % topology.atom(4))
print('All atoms: %s' % [atom for atom in topology.atoms])

print('Second residue: %s' % traj.topology.residue(1))
print('All residues: %s' % [residue for residue in traj.topology.residues])

