#!/bin/python
import mdtraj as md

#traj=md.load('../cg/sys.dcd',top='../cg/sys.pdb')
traj=md.load('../aa/sys.dcd',top='../aa/sys.pdb')

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


