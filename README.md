# cg-clusters
The scripts analyze clusters from a coarse grained trajectory frame:

python main.py -c file.dat -f frame -s file.pdb -t file.dcd

-c file.dat contains parameters (A,C and R) for a contact criterion:

dc=sqrt((A*(Ra+Rb)*0.5+D)^2)

where R is a double radius of a molecule.

OUTPUT:
1. Cluster size distribution.
2. PDB with wrapped clusters.
3. PDB with center of geometry of the clusters.
4. Radial distribution function for each cluster and each molecule type.
