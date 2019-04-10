# CG-clusters

### The repository is to determine and analyze clusters in a coarse grained simulation.

![alt text](https://github.com/Aksonik/cg-clusters/blob/master/scheme.png)


```
python main.py -c file.dat -f frame -s file.pdb -t file.dcd
```


-c file.dat contains parameters (A,C and R) for a contact criterion:

dc=A*(Ra+Rb)*0.5+D

where R is a double radius of a molecule.

## OUTPUT:
1. Cluster size distribution.
2. PDB with wrapped clusters.
3. PDB with center of geometry of the clusters.
4. Radial distribution function for each cluster and each molecule type.

## What else does it need?

mdtraj [link](http://mdtraj.org)
