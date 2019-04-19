#!/bin/python
import mdtraj as md
import numpy
import os

def molecules_types(contact):

 fc=open(contact,"r")		### read contact parameters from a file

 rns=[]		### molecule types

 for line in fc:
  w=line.split()

  if(w[0]=="Rc"):
   rns.append(w[1])

 fc.close()

 return rns

def molecules_radii(contact):

 fc=open(contact,"r")		### read contact parameters from a file

 radii={}	### molecule radii

 for line in fc:
  w=line.split()

  if(w[0]=="Rc"):
   radii[w[1]]=float(w[2])

 fc.close()

 return radii

def box(traj,dirout):
 bx=traj.unitcell_lengths[0][0]		### box
 by=traj.unitcell_lengths[0][1]
 bz=traj.unitcell_lengths[0][2]

 fi=open(str(dirout)+"/box.dat","w")
 print(bx,by,bz,file=fi)
 fi.close()

 return bx,by,bz



def cluster(traj,contact):

 fc=open(contact,"r")		### read contact parameters from a file

 radii={}
 rns=[]		### molecule types

 for line in fc:
  w=line.split()

  if(w[0]=="Ac"):
   Ac=float(w[1])
  elif(w[0]=="Dc"):
   Dc=float(w[1])
  elif(w[0]=="Rc"):
   radii[w[1]]=float(w[2])
   rns.append(w[1])

 fc.close()

 rcm={}		### matrix for number of contacts

 for rn1 in range(0,len(rns)):
  for rn2 in range(rn1,len(rns)):
   rcm[rns[rn1],rns[rn2]]=0

 c=[]		### clustered proteins
 cxyz=[]	### coordinates of clustered proteins

 bx=traj.unitcell_lengths[0][0]		### box
 by=traj.unitcell_lengths[0][1]
 bz=traj.unitcell_lengths[0][2]

 pr=[]					### remaining proteins to be clustered
 for p in range(0,traj.n_residues):
  pr.append(p)

 pl=len(pr)

### print:
### 1) protein that initiates a cluster (p0)
### 2) protein already in a cluster (p1)
### 3) protein that candidates to a cluster (p2)
### 4) clusters

### print("%3s %3s %3s %s" % ("x","x","x",c))

 cn=-1

 while pr!=[]:	### protein that initiates a cluster (one for a cluster)
   p0=pr[0]

   cn=cn+1	### cluster number

   x0=traj.xyz[0,p0,0]
   y0=traj.xyz[0,p0,1]
   z0=traj.xyz[0,p0,2]

   c.append([p0])
   pr.remove(p0)
   cxyz.append([[x0,y0,z0]])

###   print("%3i %3i %3s %s" % (p0,0,"x",c))

   pl=pl-1
   print("proteins left:",pl)

   for p1 in c[cn]:	### protein already in a cluster

     p1ndx=c[cn].index(p1)

     x1=cxyz[cn][p1ndx][0]
     y1=cxyz[cn][p1ndx][1]
     z1=cxyz[cn][p1ndx][2]

     rn1=traj.topology.atom(p1).name

     pr2rm=[]
     for p2 in pr:			### protein that candidates to a cluster

       x2=traj.xyz[0,p2,0]
       y2=traj.xyz[0,p2,1]
       z2=traj.xyz[0,p2,2]

       rn2=traj.topology.atom(p2).name

       x2=x2-bx*int((x2-x1)/bx)		### deals with unwraped trajectories
       y2=y2-by*int((y2-y1)/by)
       z2=z2-bz*int((z2-z1)/bz)

       if((x2-x1)>(0.5*bx)):		### find the closest image of the two
        x2=x2-bx
       elif((x2-x1)<(-0.5*bx)):
        x2=x2+bx

       if((y2-y1)>(0.5*by)):
        y2=y2-by
       elif((y2-y1)<(-0.5*by)):
        y2=y2+by

       if((z2-z1)>(0.5*bz)):
        z2=z2-bz
       elif((z2-z1)<(-0.5*bz)):
        z2=z2+bz

#       dxyz=numpy.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)	### closest distance
       dxyz2=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2

#       rc=numpy.sqrt((Ac*(radii[rn1]+radii[rn2])*0.5+Dc)**2)	### contact distance criterion
       rc2=(Ac*(radii[rn1]+radii[rn2])*0.5+Dc)**2

#       print("%2i-%2i (%3s-%3s) %12.6f %12.6f" % (p1,p2,rn1,rn2,rc,dxyz))

###       print("%3i %3i %3i %s" % (p0,p1,p2,c))

       if(dxyz2<rc2):		### are or are not in contact
        c[cn].append(p2)
        pr2rm.append(p2)
        cxyz[cn].append([x2,y2,z2])

        pl=pl-1
        print("proteins left:",pl)

        if (rn1,rn2) in rcm:
         rcm[rn1,rn2]=rcm[rn1,rn2]+1	### sum contacts
        else:
         rcm[rn2,rn1]=rcm[rn2,rn1]+1	### sum contacts

     for r in pr2rm:
      pr.remove(r)
   
###       print("%3i %3i %3i %s" % (p0,p1,p2,c))

 cs=sorted(c,key=len,reverse=True)		### sorted: largest -> smallest
 csxyz=sorted(cxyz,key=len,reverse=True)

 return(cs,csxyz,rcm)

def cluster_write(cs,dirout):

 f=open(str(dirout)+"/cluster.dat","w")

 for c in cs:
  print("%s:" % len(c),end=" ",file=f)
  for i in c:
   print("%i" % i,end=" ",file=f)
   if(i==(c[len(c)-1])):
    print("",file=f)		### new line

 f.close()


def contacts_write(contacts,dirout,contact):

 fc=open(contact,"r")		### read contact parameters from a file

 rns=[]		### molecule types

 for line in fc:
  w=line.split()

  if(w[0]=="Rc"):
   rns.append(w[1])

 fc.close()

 allcon=sum(contacts.values())

 f=open(str(dirout)+"/contacts.dat","w")

 for rn1 in range(0,len(rns)):
  for rn2 in range(rn1,len(rns)):
   print("%3s %3s %12i %6.3f" % (rns[rn1],rns[rn2],contacts[rns[rn1],rns[rn2]],contacts[rns[rn1],rns[rn2]]/allcon),file=f)

 f.close()



def contacts_avg(d):

 r1=[]		### type of molecule A
 r2=[]		### type of molecule B

 c=[]		### number of contacts
 p=[]		### percentage of number of contacts

 for f in d:

  fi=open("f"+str(f)+"/contacts.dat","r")

  cn=-1
  for line in fi:
   cn=cn+1

   w=line.split()

   if(f==1):
    r1.append(w[0])
    r2.append(w[1])

    c.append([float(w[2])])
    p.append([float(w[3])])
   else:
    c[cn].append(float(w[2]))
    p[cn].append(float(w[3]))

 d="contacts"
 if not os.path.exists(d):
  os.mkdir(d)

 fo=open(d+"/contacts.dat","w")

 for i in range(0,len(c)):
  print("%5s %5s %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f"
% (r1[i],r2[i],
   numpy.mean(c[i]),numpy.std(c[i]),numpy.std(c[i])/numpy.sqrt(len(c[i])),
   numpy.mean(p[i]),numpy.std(p[i]),numpy.std(p[i])/numpy.sqrt(len(p[i]))),
   file=fo)

 fo.close()
