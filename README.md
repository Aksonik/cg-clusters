# CG-clusters

### The repository is to determine and analyze clusters in a coarse grained simulation.

![alt text](https://github.com/Aksonik/cg-clusters/blob/master/scheme.png)


#### Command

```
python main.py -c file.dat -f frame -s file.pdb -t file.dcd
```

*-c file.dat* is a parameter file for the contact criterion
*-f frame* is frame number
*-s file.pdb* is a structure file (PDB)
*-t file.dcd* is a trajectory file (DCD)

#### Contact criterion

Two molecules, *A* and *B*, are assumed to be in contact if the distance between them is less than *d*:

*d* = *A<sub>c</sub>* ( *R<sub>A</sub>* + *R<sub>B</sub>* ) 0.5 + *D<sub>c</sub>*

where *A<sub>c</sub>* and *D<sub>c</sub>* are adjustable parameters, and *R* is a diameter of a molecule.

Example parameter file for the contact criterion:

```
Ac 1.0
Dc 0.7
CGA 1.5
CGB 2.0
```

where *D* and *R* values are in nanometers, CGA and CGB are atom names in the structure file.

#### Output

1. Cluster size distribution (*csd.dat*).
2. PDB with wrapped clusters (*genpdb.pdb*).
3. PDB with center of geometry of the clusters (*cog.dat*).
4. Radial distribution function for each cluster and each molecule type (*rdf*).

#### What else does it need?

MDTraj [link](http://mdtraj.org)
