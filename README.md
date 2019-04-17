# CG-clusters

### The repository is to determine and analyze clusters in coarse grained simulations.

![alt text](https://github.com/Aksonik/cg-clusters/blob/master/scheme.png)


#### Command

```
python main.py -c file -f frame -s file -t file
```

*-c file* -- a parameter file for the contact criterion

*-f frame* -- a frame number (INT)

*-s file* -- a structure file (PDB)

*-t file* -- a trajectory file (DCD)

*-fn file* -- a list of frame numbers (int)

*-bs size* -- bin size for RDF analysis (float)

#### Contact criterion

Two molecules, *A* and *B*, are assumed to be in contact if the distance between them is less than *d*:

*d* = *A<sub>c</sub>* ( *R<sub>c</sub><sup>A</sup>* + *R<sub>c</sub><sup>B</sup>* ) 0.5 + *D<sub>c</sub>*

where *A<sub>c</sub>* and *D<sub>c</sub>* are adjustable parameters, and *R<sub>c</sub>* is a diameter of a molecule.

An example of a parameter file for the contact criterion:

```
Ac 1.0
Dc 0.7
Rc CGA 1.5
Rc CGB 2.0
Rc CGC 2.5
```

where *D<sub>c</sub>* and *R<sub>c</sub>* values are in nanometers, 
*CGA* and *CGB* are atom names in the structure file.

#### Output

<ol>
<li>Cluster size distribution (<i>csd.dat, csd.png</i>).
</li>

![alt text](https://github.com/Aksonik/cg-clusters/blob/master/csd.png)

<li>PDB with wrapped clusters (<i>genpdb.pdb</i>).

Molecules belonging to different clusters have different values in the occupancy column.
The valuse are from 0.00 to 1.00, such that they can be used to color the clusters with <i>VMD</i>.
</li>

<li>PDB with center of geometry of the clusters (<i>cog.dat</i>).</li>
<li>Radial distribution function (RDF) for each cluster (<i>rdf</i>).

If not specified by the *-bs* option the bin size is equal to 0.5 nm.
</li>

<li>Number and percentage of contacts betweed different types of molecules (<i>contacts.dat</i>).</li>
</ol>


![alt text](https://github.com/Aksonik/cg-clusters/blob/master/contacts.png){width="50%"}

<!-- {:height="700px" width="400px"} -->

#### What else does it need?

MDTraj ([link](http://mdtraj.org))
