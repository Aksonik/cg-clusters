# CG-clusters

#### The repository is to determine and analyze clusters of molecules in a coarse grained simulation.

<figure>
<img src="https://github.com/Aksonik/cg-clusters/blob/master/scheme.png" width="75%" alt="" >
</figure>

<b>Figure:</b> A system of 3,000 molecules. *Left:* the molecules are colored by the size (green > red > blue). *Right:* the molecules are colored by the clusters they belong to. 

#### How does it work?

<figure>
<img src="https://github.com/Aksonik/cg-clusters/blob/master/anima.gif" width="50%" alt="" >
</figure>

<b>Figure:</b> The clusters are determined based on a contact criterion.

#### Command

```
python main.py -c file -f frame -s file -t file
```

*-c file* -- a parameter file for the contact criterion

*-f frame* -- a frame number (int)

*-fn file* -- a list of frame numbers (int)

*-s file* -- a structure file (PDB)

*-t file* -- a trajectory file (DCD)

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
*CGA*, *CGB*, and *CGC* are atom names in the structure file.

#### Output

<ol>

<li>Cluster size distribution (<i>csd.dat, csd.png</i>).</li>

![alt text](https://github.com/Aksonik/cg-clusters/blob/master/csd.png)

<li>PDB with wrapped clusters (<i>genpdb.pdb</i>).

Molecules belonging to different clusters have different values in the occupancy column.
The valuse are from 0.00 to 1.00, such that they can be used to color the clusters with <i>VMD</i>.
</li>

<li>PDB with center of geometry of the clusters (<i>cog.dat</i>).</li>
<li>Radial distribution function for each cluster (<i>rdf</i>).

<figure>
<img src="https://github.com/Aksonik/cg-clusters/blob/master/rdf.png" width="75%" alt="" >
</figure>


If not specified by the *-bs* option the bin size is equal to 0.5 nm.

An example of a file name: rdf_6_182_58_CGB.dat - a radial distribution function of 58 CGB molecules in a cluster of total 182 molecules and ID 6.
</li>

<li>Number and percentage of contacts betweed different types of molecules (<i>contacts.dat, contacts.png</i>).</li>

<figure>
<img src="https://github.com/Aksonik/cg-clusters/blob/master/contacts.png" width="75%" alt="" >
</figure>

</ol>

#### Remarks

<ol>

<li> Multiple trajectory frames are analyzed separately. Then, average, standard deviation and uncertainty is calculated.</li>

<li> If a cluster is infinite, i.e. interacts with its own images through the periodic boundary conditions, its shape is arbitrary.</li>

</ol>

#### What else does it need?

[MDTraj](http://mdtraj.org)
