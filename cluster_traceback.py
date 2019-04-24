#!/bin/python

import sys
import numpy

def cluster_traceback():

### a cluster that is going to be traced back

 ff=6000
 ll=0

 fi=open("f"+str(ff)+"/cluster.dat","r")

 l=-1
 for line in fi:
  l=l+1
  if(l==0):
   w=line.split(" ")
   c=w[1:len(w)-1]	### not the first (size) and not the last (new line)
   cc=c
   break

 fi.close()



### compare with clusters from the previous frames

 for f in range(6000,0,-50):
  fi=open("f"+str(f)+"/cluster.dat","r")
  fo=open("f"+str(f)+"/cluster_traceback.dat","w")

  permax=0

  for line in fi:

   w=line.split(" ")
   c=w[1:len(w)-1]	### not the first (size) and not the last (new line)

   co=(set(c)&set(cc))	### molecules common for the previous and current clusters
   per=float(len(co))/float(len(cc))

   if(per>permax):
    permax=per
    cmax=c
    ccmax=c

#  print("previous cluster:",cc)
#  print("current cluster:",ccmax)
#  print("common molecules: %12.6f" % permax)

  cc=ccmax

  print("%9i %6.3f %s" % (f,permax,cc),file=fo)	### this format needs to be fixed (only numbers)

  ### for now:
  ### sed -i '{s:\[::g}' f*/cluster_traceback.dat 
  ### sed -i '{s:\]::g}' f*/cluster_traceback.dat 
  ### sed -i '{s:,::g}' f*/cluster_traceback.dat 
  ### sed -i "s/'//g" f*/cluster_traceback.dat 

  fi.close()
  fo.close()

def cluster_traceback_write(n):

### read molecules in the cluster

 cluster=[]
 fi=open("f"+str(n)+"/cluster_traceback.dat","r")
 for line in fi:
  w=line.split()
  c=w[2:]
  for j in c:
   cluster.append(int(j)+1)
 fi.close()

### recognize them in the pdb file

 fi=open("f"+str(n)+"/genpdb.pdb","r")
 fo=open("f"+str(n)+"/cluster_traceback.pdb","w")
 fovmd=open("f"+str(n)+"/cluster_traceback_vmd.pdb","w")

 for line in fi:
  w=line.split()
  atomtype=str(line[0:6].strip())
  if(atomtype=="CRYST1"):
   print(line,file=fo)
   print(line,file=fovmd)
  if(atomtype=="ATOM"):

   atomnumber=int(line[6:11].strip())
   atomname=str(line[12:16].strip())
   alterloc=str(line[16:17].strip())
   resname=str(line[17:20].strip())
   chain=str(line[21:22].strip())
   resseqnum=int(line[22:26].strip())
   resinsertion=str(line[26:27].strip())
   x=float(line[30:38].strip())
   y=float(line[38:46].strip())
   z=float(line[46:54].strip())
   occupancy=float(line[54:60].strip())
   tempfac=float(line[60:66].strip())
#   element=str(line[76:78].strip())
#   charge=str(line[78:80].strip())
   segname=str(line[72:76].strip())

### vmd doesn't work for trajectories with varying number of atoms.
### it also doesn't work for varying occupancy - reads it only once.

   if atomnumber not in cluster:
    x=0.000
    y=0.000
    z=0.000
   else:
    print("%-6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s" % 
(atomtype,atomnumber,atomname,alterloc,resname,chain,resseqnum,resinsertion,x,y,z,
occupancy,tempfac,segname),file=fo)

   print("%-6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s" % 
(atomtype,atomnumber,atomname,alterloc,resname,chain,resseqnum,resinsertion,x,y,z,
occupancy,tempfac,segname),file=fovmd)

 fi.close()
 fo.close()
 fovmd.close()
