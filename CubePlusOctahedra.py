# This code defines a superstrcuture of Cubic and Octahedral Colloidal Clusters
# Author: Zanjani research group: http://cec.miamioh.edu/zanjani-research 
# Details about the structure are presented in the following reference:
# Aryana, K., Stahley, J. B., Parvez, N., Kim, K. and Zanjani, M. B. (2019) Superstructures of multielement colloidal molecules: efficient pathways to construct reconfigurable photonic and phononic crystals. Adv. Theory Simul. p. 1800198.
# To visualize the configuration, you can load the output text file into the Ovito (or similar software), and perform an Affine Transformation that scale everything with a factor of 1e5
# Following the previous step in Ovito, you can select the following radii for various types
# raduis of type 1 for visualization in ovito = 2
# raduis of type 2 for visualization in ovito = 3.75
# raduis of type 3 for visualization in ovito = 2
# raduis of type 4 for visualization in ovito = 1.46

import random
import numpy as np
import math

#particle sizes are chosen for one specific scenario, they can be freely changed according to geometric constraints discussed in the article

dB=400e-7 #defines the size of the vertex particles of a cubic cluster in CGS units
dBeq=dB+11e-7 # offset for DNA-mediated interactions
poldensity=1.05 # core density in g/cm^3, example: polystyrine
massB=poldensity*math.pi*(dB**3)/6.0
b=dBeq
rEB=1.0
dE=dB #size of vertex particle of octahedral units 
massE=poldensity*math.pi*(dE**3)/6.0
a=dBeq*1.91 #lattice parameter calculated based on geometric constraints explained in the article

dA=1.875*dB # size of the center particle of the octahedral unit
massA=poldensity*math.pi*(dA**3)/6.0
e=dBeq
dD=0.732 * dB #size of the center particle of the cubic unit
massD=poldensity*math.pi*(dD**3)/6

atompercell=4*9+4*7
xs=[0.00] * atompercell
ys=[0.00] * atompercell
zs=[0.00] * atompercell
types=[0] * atompercell
molid0=[0] * atompercell

#define (4) octahedra
inde=-1
xcc=a
ycc=0.0
zcc=0.0

xs[inde+7]=xcc
ys[inde+7]=ycc
zs[inde+7]=zcc
types[inde+7]=2
molid0[inde+7]=1

xs[inde+1]=xcc-e/math.sqrt(2.0)
ys[inde+1]=ycc
zs[inde+1]=zcc
types[inde+1]=1
molid0[inde+1]=1

xs[inde+2]=xcc
ys[inde+2]=ycc+e/math.sqrt(2.0)
zs[inde+2]=zcc
types[inde+2]=1
molid0[inde+2]=1

xs[inde+3]=xcc+e/math.sqrt(2.0)
ys[inde+3]=ycc
zs[inde+3]=zcc
types[inde+3]=1
molid0[inde+3]=1

xs[inde+4]=xcc
ys[inde+4]=ycc-e/math.sqrt(2.0)
zs[inde+4]=zcc
types[inde+4]=1
molid0[inde+4]=1

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc+e/math.sqrt(2.0)
types[inde+5]=1
molid0[inde+5]=1

xs[inde+6]=xcc
ys[inde+6]=ycc
zs[inde+6]=zcc-e/math.sqrt(2.0)
types[inde+6]=1
molid0[inde+6]=1

inde=6
xcc=0.0
ycc=a
zcc=0.0

xs[inde+7]=xcc
ys[inde+7]=ycc
zs[inde+7]=zcc
types[inde+7]=2
molid0[inde+7]=2

xs[inde+1]=xcc-e/math.sqrt(2.0)
ys[inde+1]=ycc
zs[inde+1]=zcc
types[inde+1]=1
molid0[inde+1]=2

xs[inde+2]=xcc
ys[inde+2]=ycc+e/math.sqrt(2.0)
zs[inde+2]=zcc
types[inde+2]=1
molid0[inde+2]=2

xs[inde+3]=xcc+e/math.sqrt(2.0)
ys[inde+3]=ycc
zs[inde+3]=zcc
types[inde+3]=1
molid0[inde+3]=2

xs[inde+4]=xcc
ys[inde+4]=ycc-e/math.sqrt(2.0)
zs[inde+4]=zcc
types[inde+4]=1
molid0[inde+4]=2

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc+e/math.sqrt(2.0)
types[inde+5]=1
molid0[inde+5]=2

xs[inde+6]=xcc
ys[inde+6]=ycc
zs[inde+6]=zcc-e/math.sqrt(2.0)
types[inde+6]=1
molid0[inde+6]=2

inde=13
xcc=0.0
ycc=0.0
zcc=a

xs[inde+7]=xcc
ys[inde+7]=ycc
zs[inde+7]=zcc
types[inde+7]=2
molid0[inde+7]=3

xs[inde+1]=xcc-e/math.sqrt(2.0)
ys[inde+1]=ycc
zs[inde+1]=zcc
types[inde+1]=1
molid0[inde+1]=3

xs[inde+2]=xcc
ys[inde+2]=ycc+e/math.sqrt(2.0)
zs[inde+2]=zcc
types[inde+2]=1
molid0[inde+2]=3

xs[inde+3]=xcc+e/math.sqrt(2.0)
ys[inde+3]=ycc
zs[inde+3]=zcc
types[inde+3]=1
molid0[inde+3]=3

xs[inde+4]=xcc
ys[inde+4]=ycc-e/math.sqrt(2.0)
zs[inde+4]=zcc
types[inde+4]=1
molid0[inde+4]=3

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc+e/math.sqrt(2.0)
types[inde+5]=1
molid0[inde+5]=3

xs[inde+6]=xcc
ys[inde+6]=ycc
zs[inde+6]=zcc-e/math.sqrt(2.0)
types[inde+6]=1
molid0[inde+6]=3

inde=20
xcc=a
ycc=a
zcc=a

xs[inde+7]=xcc
ys[inde+7]=ycc
zs[inde+7]=zcc
types[inde+7]=2
molid0[inde+7]=4

xs[inde+1]=xcc-e/math.sqrt(2.0)
ys[inde+1]=ycc
zs[inde+1]=zcc
types[inde+1]=1
molid0[inde+1]=4

xs[inde+2]=xcc
ys[inde+2]=ycc+e/math.sqrt(2.0)
zs[inde+2]=zcc
types[inde+2]=1
molid0[inde+2]=4

xs[inde+3]=xcc+e/math.sqrt(2.0)
ys[inde+3]=ycc
zs[inde+3]=zcc
types[inde+3]=1
molid0[inde+3]=4

xs[inde+4]=xcc
ys[inde+4]=ycc-e/math.sqrt(2.0)
zs[inde+4]=zcc
types[inde+4]=1
molid0[inde+4]=4

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc+e/math.sqrt(2.0)
types[inde+5]=1
molid0[inde+5]=4

xs[inde+6]=xcc
ys[inde+6]=ycc
zs[inde+6]=zcc-e/math.sqrt(2.0)
types[inde+6]=1
molid0[inde+6]=4

#define (4) Cubes
inde=27
xcc=0.0
ycc=0.0
zcc=0.0

xs[inde+9]=xcc
ys[inde+9]=ycc
zs[inde+9]=zcc
types[inde+9]=4
molid0[inde+9]=5
xs[inde+1]=xcc-b/2.0
ys[inde+1]=ycc-b/2.0
zs[inde+1]=zcc-b/2.0
types[inde+1]=3
molid0[inde+1]=5
xs[inde+2]=xcc+b/2.0
ys[inde+2]=ycc-b/2.0
zs[inde+2]=zcc-b/2.0
types[inde+2]=3
molid0[inde+2]=5
xs[inde+3]=xcc+b/2.0
ys[inde+3]=ycc+b/2.0
zs[inde+3]=zcc-b/2.0
types[inde+3]=3
molid0[inde+3]=5
xs[inde+4]=xcc-b/2.0
ys[inde+4]=ycc+b/2.0
zs[inde+4]=zcc-b/2.0
types[inde+4]=3
molid0[inde+4]=5

xs[inde+5]=xcc-b/2.0
ys[inde+5]=ycc-b/2.0
zs[inde+5]=zcc+b/2.0
types[inde+5]=3
molid0[inde+5]=5
xs[inde+6]=xcc+b/2.0
ys[inde+6]=ycc-b/2.0
zs[inde+6]=zcc+b/2.0
types[inde+6]=3
molid0[inde+6]=5
xs[inde+7]=xcc+b/2.0
ys[inde+7]=ycc+b/2.0
zs[inde+7]=zcc+b/2.0
types[inde+7]=3
molid0[inde+7]=5
xs[inde+8]=xcc-b/2.0
ys[inde+8]=ycc+b/2.0
zs[inde+8]=zcc+b/2.0
types[inde+8]=3
molid0[inde+8]=5

inde=36
xcc=a
ycc=a
zcc=0.0

xs[inde+9]=xcc
ys[inde+9]=ycc
zs[inde+9]=zcc
types[inde+9]=4
molid0[inde+9]=6
xs[inde+1]=xcc-b/2.0
ys[inde+1]=ycc-b/2.0
zs[inde+1]=zcc-b/2.0
types[inde+1]=3
molid0[inde+1]=6
xs[inde+2]=xcc+b/2.0
ys[inde+2]=ycc-b/2.0
zs[inde+2]=zcc-b/2.0
types[inde+2]=3
molid0[inde+2]=6
xs[inde+3]=xcc+b/2.0
ys[inde+3]=ycc+b/2.0
zs[inde+3]=zcc-b/2.0
types[inde+3]=3
molid0[inde+3]=6
xs[inde+4]=xcc-b/2.0
ys[inde+4]=ycc+b/2.0
zs[inde+4]=zcc-b/2.0
types[inde+4]=3
molid0[inde+4]=6

xs[inde+5]=xcc-b/2.0
ys[inde+5]=ycc-b/2.0
zs[inde+5]=zcc+b/2.0
types[inde+5]=3
molid0[inde+5]=6
xs[inde+6]=xcc+b/2.0
ys[inde+6]=ycc-b/2.0
zs[inde+6]=zcc+b/2.0
types[inde+6]=3
molid0[inde+6]=6
xs[inde+7]=xcc+b/2.0
ys[inde+7]=ycc+b/2.0
zs[inde+7]=zcc+b/2.0
types[inde+7]=3
molid0[inde+7]=6
xs[inde+8]=xcc-b/2.0
ys[inde+8]=ycc+b/2.0
zs[inde+8]=zcc+b/2.0
types[inde+8]=3
molid0[inde+8]=6

inde=45
xcc=a
ycc=0.0
zcc=a

xs[inde+9]=xcc
ys[inde+9]=ycc
zs[inde+9]=zcc
types[inde+9]=4
molid0[inde+9]=7
xs[inde+1]=xcc-b/2.0
ys[inde+1]=ycc-b/2.0
zs[inde+1]=zcc-b/2.0
types[inde+1]=3
molid0[inde+1]=7
xs[inde+2]=xcc+b/2.0
ys[inde+2]=ycc-b/2.0
zs[inde+2]=zcc-b/2.0
types[inde+2]=3
molid0[inde+2]=7
xs[inde+3]=xcc+b/2.0
ys[inde+3]=ycc+b/2.0
zs[inde+3]=zcc-b/2.0
types[inde+3]=3
molid0[inde+3]=7
xs[inde+4]=xcc-b/2.0
ys[inde+4]=ycc+b/2.0
zs[inde+4]=zcc-b/2.0
types[inde+4]=3
molid0[inde+4]=7

xs[inde+5]=xcc-b/2.0
ys[inde+5]=ycc-b/2.0
zs[inde+5]=zcc+b/2.0
types[inde+5]=3
molid0[inde+5]=7
xs[inde+6]=xcc+b/2.0
ys[inde+6]=ycc-b/2.0
zs[inde+6]=zcc+b/2.0
types[inde+6]=3
molid0[inde+6]=7
xs[inde+7]=xcc+b/2.0
ys[inde+7]=ycc+b/2.0
zs[inde+7]=zcc+b/2.0
types[inde+7]=3
molid0[inde+7]=7
xs[inde+8]=xcc-b/2.0
ys[inde+8]=ycc+b/2.0
zs[inde+8]=zcc+b/2.0
types[inde+8]=3
molid0[inde+8]=7

inde=54
xcc=0.0
ycc=a
zcc=a

xs[inde+9]=xcc
ys[inde+9]=ycc
zs[inde+9]=zcc
types[inde+9]=4
molid0[inde+9]=8
xs[inde+1]=xcc-b/2.0
ys[inde+1]=ycc-b/2.0
zs[inde+1]=zcc-b/2.0
types[inde+1]=3
molid0[inde+1]=8
xs[inde+2]=xcc+b/2.0
ys[inde+2]=ycc-b/2.0
zs[inde+2]=zcc-b/2.0
types[inde+2]=3
molid0[inde+2]=8
xs[inde+3]=xcc+b/2.0
ys[inde+3]=ycc+b/2.0
zs[inde+3]=zcc-b/2.0
types[inde+3]=3
molid0[inde+3]=8
xs[inde+4]=xcc-b/2.0
ys[inde+4]=ycc+b/2.0
zs[inde+4]=zcc-b/2.0
types[inde+4]=3
molid0[inde+4]=8

xs[inde+5]=xcc-b/2.0
ys[inde+5]=ycc-b/2.0
zs[inde+5]=zcc+b/2.0
types[inde+5]=3
molid0[inde+5]=8
xs[inde+6]=xcc+b/2.0
ys[inde+6]=ycc-b/2.0
zs[inde+6]=zcc+b/2.0
types[inde+6]=3
molid0[inde+6]=8
xs[inde+7]=xcc+b/2.0
ys[inde+7]=ycc+b/2.0
zs[inde+7]=zcc+b/2.0
types[inde+7]=3
molid0[inde+7]=8
xs[inde+8]=xcc-b/2.0
ys[inde+8]=ycc+b/2.0
zs[inde+8]=zcc+b/2.0
types[inde+8]=3
molid0[inde+8]=8


maxx=3 #number of replica of unit cell in x direction
maxy=3 #number of replica of unit cell in y direction
maxz=3 #number of replica of unit cell in z direction
nt=0
xcenter=0.0
ycenter=0.0
zcenter=0.0
cellcount=0
molpercell=8
bondpercellAE=4*18
bondpercellB=4*20
bondpercell=bondpercellAE+bondpercellB

totalatoms= atompercell * maxx * maxy * maxz
xpos = [0.00] * totalatoms
ypos = [0.00] * totalatoms
zpos = [0.00] * totalatoms
typepos = [0] * totalatoms
molidcore= [0] * totalatoms

tp1count=0
tp3count=0
idoftp1=[]
idoftp3=[]
for nx in range(maxx):
  for ny in range(maxy):
    for nz in range(maxz):
      cellcount += 1
      for i in range(atompercell):
        xpos[nt]=xs[i]+nx*a*2
        ypos[nt]=ys[i]+ny*a*2
        zpos[nt]=zs[i]+nz*a*2
        typepos[nt]=types[i]
        if typepos[nt]==1:
          tp1count += 1
          idoftp1.append(nt)
        if typepos[nt]==3:
          tp3count += 1
          idoftp3.append(nt)
        xcenter=xcenter+xpos[nt]
        ycenter=ycenter+ypos[nt]
        zcenter=zcenter+zpos[nt]
        molidcore[nt]=(cellcount-1)*molpercell+molid0[i]
        nt += 1

molseed=cellcount*molpercell
xcenter=xcenter/float(nt)
ycenter=ycenter/float(nt)
zcenter=zcenter/float(nt)

for i in range(nt):
   xpos[i] -= xcenter
   ypos[i] -= ycenter
   zpos[i] -= zcenter

corexmin= min(xpos)
corexmax= max(xpos)

xlo=corexmin
xhi=xlo+(maxx)*2*a
zoff=0.0
boxx=xhi-xlo
boxy=boxx
boxz=boxx

bondnum=bondpercell*cellcount
atomtypes=4
bondtypescount=4

flammps = open("CUBEOCT_initLAMMPS_3by3.txt", "w") #lammps data file preparation
flammps.write("LAMMPS input data for multiple clusters\r\n\r\n")
flammps.write("%d atoms\r\n" % nt)
flammps.write("%d bonds\r\n" % bondnum)

flammps.write("%d atom types\r\n" % atomtypes)
flammps.write("%d bond types\r\n" % bondtypescount)

flammps.write("%f %f xlo xhi\r\n" % (xlo, xhi))
flammps.write("%f %f ylo yhi\r\n" % (xlo, xhi))
flammps.write("%f %f zlo zhi\r\n" % (xlo, xhi))

flammps.write("\nMasses\n\n")
flammps.write("%d %fe-14 \r\n" % (1, massE*1e14))
flammps.write("%d %fe-14 \r\n" % (2, massA*1e14))
flammps.write("%d %fe-14 \r\n" % (3, massB*1e14))
flammps.write("%d %fe-14 \r\n" % (4, massD*1e14))


#Atoms
flammps.write("\nAtoms\r\n\n")
for i0 in range(nt):
  flammps.write("%d %d %d %f %f %f\r\n" % (i0+1, molidcore[i0], typepos[i0], xpos[i0], ypos[i0], zpos[i0]))
#Bonds
flammps.write("\nBonds\r\n\n")
bn=0
for i in range(cellcount):
  bn = i * bondpercell
  molstart= i * atompercell
  for j in range(4): #for octahedra
    bondtype=1
    at_n1=molstart+j*7+1
    at_n2=molstart+j*7+2
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+1, bondtype,at_n1,at_n2))

    at_n1=molstart+j*7+2
    at_n2=molstart+j*7+3
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+2, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+3
    at_n2=molstart+j*7+4
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+3, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+1
    at_n2=molstart+j*7+4
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+4, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+1
    at_n2=molstart+j*7+5
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+5, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+2
    at_n2=molstart+j*7+5
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+6, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+3
    at_n2=molstart+j*7+5
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+7, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+4
    at_n2=molstart+j*7+5
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+8, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+1
    at_n2=molstart+j*7+6
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+9, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+2
    at_n2=molstart+j*7+6
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+10, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+3
    at_n2=molstart+j*7+6
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+11, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+4
    at_n2=molstart+j*7+6
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+12, bondtype,at_n1,at_n2))
    
    bondtype=2
    at_n1=molstart+j*7+1
    at_n2=molstart+j*7+7
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+13, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+2
    at_n2=molstart+j*7+7
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+14, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+3
    at_n2=molstart+j*7+7
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+15, bondtype,at_n1,at_n2))
  
    at_n1=molstart+j*7+4
    at_n2=molstart+j*7+7
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+16, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+5
    at_n2=molstart+j*7+7
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+17, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*7+6
    at_n2=molstart+j*7+7
    flammps.write("%d %d %d %d\r\n" % (bn+j*18+18, bondtype,at_n1,at_n2))
    
  #cubes
  for j in range(4): 
    bondtype=1
    at_n1=molstart+4*7+j*9+1
    at_n2=molstart+4*7+j*9+2
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+1, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+2
    at_n2=molstart+4*7+j*9+3
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+2, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+3
    at_n2=molstart+4*7+j*9+4
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+3, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+1
    at_n2=molstart+4*7+j*9+4
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+4, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+5
    at_n2=molstart+4*7+j*9+6
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+5, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+6
    at_n2=molstart+4*7+j*9+7
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+6, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+7
    at_n2=molstart+4*7+j*9+8
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+7, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+5
    at_n2=molstart+4*7+j*9+8
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+8, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+1
    at_n2=molstart+4*7+j*9+5
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+9, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+2
    at_n2=molstart+4*7+j*9+6
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+10, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+3
    at_n2=molstart+4*7+j*9+7
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+11, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+4
    at_n2=molstart+4*7+j*9+8
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+12, bondtype,at_n1,at_n2))
    
    bondtype=4
    at_n1=molstart+4*7+j*9+1
    at_n2=molstart+4*7+j*9+9
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+13, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+2
    at_n2=molstart+4*7+j*9+9
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+14, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+3
    at_n2=molstart+4*7+j*9+9
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+15, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+4
    at_n2=molstart+4*7+j*9+9
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+16, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+5
    at_n2=molstart+4*7+j*9+9
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+17, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+6
    at_n2=molstart+4*7+j*9+9
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+18, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+7
    at_n2=molstart+4*7+j*9+9
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+19, bondtype,at_n1,at_n2))
    
    at_n1=molstart+4*7+j*9+8
    at_n2=molstart+4*7+j*9+9
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*20+20, bondtype,at_n1,at_n2))
    
flammps.close()
print("done with Lammps read_data file")
