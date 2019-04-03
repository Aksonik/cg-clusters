#!/bin/python

def cog(clust,cxyz,traj):

 f=open("cog.pdb","w")

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


  r="COG"
  n=1
  print("%6s%5i%5s%4s%6i    %8.3f%8.3f%8.3f%s%s" % 
   ("ATOM  ",n,r,r,n,cogx/len(c)*10.0,cogy/len(c)*10.0,cogz/len(c)*10.0,"  1.00  1.00      ","PRO0"),file=f)	### [nm] -> [A]

 f.close()

 return cog
