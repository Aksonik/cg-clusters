# CG-clusters

### The repository is to determine and analyze clusters in coarse grained simulations.

![alt text](https://github.com/Aksonik/cg-clusters/blob/master/scheme.png)


#### Command

```
python main.py -c file -f frame -s file -t file
```

*-c file.dat* is a parameter file for the contact criterion

*-f frame* is a frame number (TXT)

*-s file.pdb* is a structure file (PDB)

*-t file.dcd* is a trajectory file (DCD)

#### Contact criterion

Two molecules, *A* and *B*, are assumed to be in contact if the distance between them is less than *d*:

```
*d* = *A<sub>c</sub>* ( *R<sub>c</sub><sup>A</sup>* + *R<sub>c</sub><sup>B</sup>* ) 0.5 + *D<sub>c</sub>*
```

where *A<sub>c</sub>* and *D<sub>c</sub>* are adjustable parameters, and *R<sub>c</sub>* is a diameter of a molecule.

An example parameter file for the contact criterion:

```
Ac 1.0
Dc 0.7
Rc CGA 1.5
Rc CGB 2.0
```

where *D<sub>c</sub>* and *R<sub>c</sub>* values are in nanometers, 
*CGA* and *CGB* are atom names in the structure file.

#### Output

1. Cluster size distribution (*csd.dat*).
2. PDB with wrapped clusters (*genpdb.pdb*).

Molecules belonging to different clusters have different values in the occupancy column.
The valuse are from 0.00 to 1.00, such that they can be used to color the clusters with *VMD*.

3. PDB with center of geometry of the clusters (*cog.dat*).
4. Radial distribution function for each cluster (*rdf*).

#### What else does it need?

MDTraj [link](http://mdtraj.org)
