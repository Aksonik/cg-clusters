#!/bin/python

def cog(clust,cxyz,traj):

 cog=[]

 for c in clust:

  cogx=0
  cogy=0
  cogz=0

  for p in c:

   x=cxyz[clust.index(c)][c.index(p)][0]	### [nm]
   y=cxyz[clust.index(c)][c.index(p)][1]
   z=cxyz[clust.index(c)][c.index(p)][2]

   cogx=cogx+x
   cogy=cogy+y
   cogz=cogz+z

  cog.append([cogx/len(c),cogy/len(c),cogz/len(c)])

 return cog



def cog_write(cog,dirout):

 f=open(str(dirout)+"/cog.pdb","w")

 r="COG"
 n=1

 cn=0

 for c in cog:

  cn=cn+1
  o=(cn%10)*0.1
  b=1.00

  cogx=c[0]
  cogy=c[1]
  cogz=c[2]

  print("%6s%5i%5s%4s%6i    %8.3f%8.3f%8.3f  %4.2f  %4.2f      %s" % 
  ("ATOM  ",n,r,r,n,cogx*10.0,cogy*10.0,cogz*10.0,o,b,"PRO0"),file=f)	### [nm] -> [A]

 f.close()
