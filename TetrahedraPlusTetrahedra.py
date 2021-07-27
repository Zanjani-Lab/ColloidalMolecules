# This MATLAB code defines a superstructure of two types of Tetrahedral Colloidal Clusters
# Details about the structure are presented in the following reference:
# Aryana, K., Stahley, J. B., Parvez, N., Kim, K. and Zanjani, M. B. (2019) Superstructures of multielement colloidal molecules: efficient pathways to construct reconfigurable photonic and phononic crystals. Adv. Theory Simul. p. 1800198.
# Author: Zanjani research group: http://cec.miamioh.edu/zanjani-research 
# For more details about the code setup, visualization in Ovito, and performing MD simulations refer to the code provided for Cubic + Octahedral units
#particle sizes are chosen arbitrarily in this code for one sample scenario, they can be freely changed according to geometric constraints discussed in the article

import random
import numpy as np
import math

dB=400e-7
dBeq=dB+11e-7
poldensity=1.05
massB=poldensity*math.pi*(dB**3)/6.0


b=dBeq*1.47

a=math.sqrt(2.0)*(dBeq+b)
dA=1.0*dB
dE=1.15*dB
e=dBeq*0.84
dD=0.85*dB
dBn=dB
b=math.sqrt(2.0/3.0)*(dD+dBn)
massE=poldensity*math.pi*(dE**3)/6.0
massA=poldensity*math.pi*(dA**3)/6.0
massD=poldensity*math.pi*(dD**3)/6.0

atompercell=4*5+4*5
xs=[0.00] * atompercell
ys=[0.00] * atompercell
zs=[0.00] * atompercell
types=[0] * atompercell
molid0=[0] * atompercell


inde=-1
xcc=0.25*a
ycc=0.25*a
zcc=0.25*a

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc
types[inde+5]=2
molid0[inde+5]=1

xs[inde+1]=xcc-e/math.sqrt(8.0)
ys[inde+1]=ycc+e/math.sqrt(8.0)
zs[inde+1]=zcc+e/math.sqrt(8.0)
types[inde+1]=1
molid0[inde+1]=1

xs[inde+2]=xcc-e/math.sqrt(8.0)
ys[inde+2]=ycc-e/math.sqrt(8.0)
zs[inde+2]=zcc-e/math.sqrt(8.0)
types[inde+2]=1
molid0[inde+2]=1

xs[inde+3]=xcc+e/math.sqrt(8.0)
ys[inde+3]=ycc-e/math.sqrt(8.0)
zs[inde+3]=zcc+e/math.sqrt(8.0)
types[inde+3]=1
molid0[inde+3]=1

xs[inde+4]=xcc+e/math.sqrt(8.0)
ys[inde+4]=ycc+e/math.sqrt(8.0)
zs[inde+4]=zcc-e/math.sqrt(8.0)
types[inde+4]=1
molid0[inde+4]=1

inde=4
xcc=0.75*a
ycc=0.75*a
zcc=0.25*a

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc
types[inde+5]=2
molid0[inde+5]=2

xs[inde+1]=xcc-e/math.sqrt(8.0)
ys[inde+1]=ycc+e/math.sqrt(8.0)
zs[inde+1]=zcc+e/math.sqrt(8.0)
types[inde+1]=1
molid0[inde+1]=2

xs[inde+2]=xcc-e/math.sqrt(8.0)
ys[inde+2]=ycc-e/math.sqrt(8.0)
zs[inde+2]=zcc-e/math.sqrt(8.0)
types[inde+2]=1
molid0[inde+2]=2

xs[inde+3]=xcc+e/math.sqrt(8.0)
ys[inde+3]=ycc-e/math.sqrt(8.0)
zs[inde+3]=zcc+e/math.sqrt(8.0)
types[inde+3]=1
molid0[inde+3]=2

xs[inde+4]=xcc+e/math.sqrt(8.0)
ys[inde+4]=ycc+e/math.sqrt(8.0)
zs[inde+4]=zcc-e/math.sqrt(8.0)
types[inde+4]=1
molid0[inde+4]=2

inde=9
xcc=0.25*a
ycc=0.75*a
zcc=0.75*a

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc
types[inde+5]=2
molid0[inde+5]=3

xs[inde+1]=xcc-e/math.sqrt(8.0)
ys[inde+1]=ycc+e/math.sqrt(8.0)
zs[inde+1]=zcc+e/math.sqrt(8.0)
types[inde+1]=1
molid0[inde+1]=3

xs[inde+2]=xcc-e/math.sqrt(8.0)
ys[inde+2]=ycc-e/math.sqrt(8.0)
zs[inde+2]=zcc-e/math.sqrt(8.0)
types[inde+2]=1
molid0[inde+2]=3

xs[inde+3]=xcc+e/math.sqrt(8.0)
ys[inde+3]=ycc-e/math.sqrt(8.0)
zs[inde+3]=zcc+e/math.sqrt(8.0)
types[inde+3]=1
molid0[inde+3]=3

xs[inde+4]=xcc+e/math.sqrt(8.0)
ys[inde+4]=ycc+e/math.sqrt(8.0)
zs[inde+4]=zcc-e/math.sqrt(8.0)
types[inde+4]=1
molid0[inde+4]=3

inde=14
xcc=0.75*a
ycc=0.25*a
zcc=0.75*a

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc
types[inde+5]=2
molid0[inde+5]=4

xs[inde+1]=xcc-e/math.sqrt(8.0)
ys[inde+1]=ycc+e/math.sqrt(8.0)
zs[inde+1]=zcc+e/math.sqrt(8.0)
types[inde+1]=1
molid0[inde+1]=4

xs[inde+2]=xcc-e/math.sqrt(8.0)
ys[inde+2]=ycc-e/math.sqrt(8.0)
zs[inde+2]=zcc-e/math.sqrt(8.0)
types[inde+2]=1
molid0[inde+2]=4

xs[inde+3]=xcc+e/math.sqrt(8.0)
ys[inde+3]=ycc-e/math.sqrt(8.0)
zs[inde+3]=zcc+e/math.sqrt(8.0)
types[inde+3]=1
molid0[inde+3]=4

xs[inde+4]=xcc+e/math.sqrt(8.0)
ys[inde+4]=ycc+e/math.sqrt(8.0)
zs[inde+4]=zcc-e/math.sqrt(8.0)
types[inde+4]=1
molid0[inde+4]=4

inde=19
xcc=0.0
ycc=0.0
zcc=0.0

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc
types[inde+5]=4
molid0[inde+5]=5

xs[inde+1]=xcc-b/math.sqrt(8.0)
ys[inde+1]=ycc+b/math.sqrt(8.0)
zs[inde+1]=zcc+b/math.sqrt(8.0)
types[inde+1]=3
molid0[inde+1]=5

xs[inde+2]=xcc-b/math.sqrt(8.0)
ys[inde+2]=ycc-b/math.sqrt(8.0)
zs[inde+2]=zcc-b/math.sqrt(8.0)
types[inde+2]=3
molid0[inde+2]=5

xs[inde+3]=xcc+b/math.sqrt(8.0)
ys[inde+3]=ycc-b/math.sqrt(8.0)
zs[inde+3]=zcc+b/math.sqrt(8.0)
types[inde+3]=3
molid0[inde+3]=5

xs[inde+4]=xcc+b/math.sqrt(8.0)
ys[inde+4]=ycc+b/math.sqrt(8.0)
zs[inde+4]=zcc-b/math.sqrt(8.0)
types[inde+4]=3
molid0[inde+4]=5

inde=24
xcc=0.5*a
ycc=0.5*a
zcc=0.0

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc
types[inde+5]=4
molid0[inde+5]=6

xs[inde+1]=xcc-b/math.sqrt(8.0)
ys[inde+1]=ycc+b/math.sqrt(8.0)
zs[inde+1]=zcc+b/math.sqrt(8.0)
types[inde+1]=3
molid0[inde+1]=6

xs[inde+2]=xcc-b/math.sqrt(8.0)
ys[inde+2]=ycc-b/math.sqrt(8.0)
zs[inde+2]=zcc-b/math.sqrt(8.0)
types[inde+2]=3
molid0[inde+2]=6

xs[inde+3]=xcc+b/math.sqrt(8.0)
ys[inde+3]=ycc-b/math.sqrt(8.0)
zs[inde+3]=zcc+b/math.sqrt(8.0)
types[inde+3]=3
molid0[inde+3]=6

xs[inde+4]=xcc+b/math.sqrt(8.0)
ys[inde+4]=ycc+b/math.sqrt(8.0)
zs[inde+4]=zcc-b/math.sqrt(8.0)
types[inde+4]=3
molid0[inde+4]=6

inde=29
xcc=0.0
ycc=0.5*a
zcc=0.5*a

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc
types[inde+5]=4
molid0[inde+5]=7

xs[inde+1]=xcc-b/math.sqrt(8.0)
ys[inde+1]=ycc+b/math.sqrt(8.0)
zs[inde+1]=zcc+b/math.sqrt(8.0)
types[inde+1]=3
molid0[inde+1]=7

xs[inde+2]=xcc-b/math.sqrt(8.0)
ys[inde+2]=ycc-b/math.sqrt(8.0)
zs[inde+2]=zcc-b/math.sqrt(8.0)
types[inde+2]=3
molid0[inde+2]=7

xs[inde+3]=xcc+b/math.sqrt(8.0)
ys[inde+3]=ycc-b/math.sqrt(8.0)
zs[inde+3]=zcc+b/math.sqrt(8.0)
types[inde+3]=3
molid0[inde+3]=7

xs[inde+4]=xcc+b/math.sqrt(8.0)
ys[inde+4]=ycc+b/math.sqrt(8.0)
zs[inde+4]=zcc-b/math.sqrt(8.0)
types[inde+4]=3
molid0[inde+4]=7

inde=34
xcc=0.5*a
ycc=0.0
zcc=0.5*a

xs[inde+5]=xcc
ys[inde+5]=ycc
zs[inde+5]=zcc
types[inde+5]=4
molid0[inde+5]=8

xs[inde+1]=xcc-b/math.sqrt(8.0)
ys[inde+1]=ycc+b/math.sqrt(8.0)
zs[inde+1]=zcc+b/math.sqrt(8.0)
types[inde+1]=3
molid0[inde+1]=8

xs[inde+2]=xcc-b/math.sqrt(8.0)
ys[inde+2]=ycc-b/math.sqrt(8.0)
zs[inde+2]=zcc-b/math.sqrt(8.0)
types[inde+2]=3
molid0[inde+2]=8

xs[inde+3]=xcc+b/math.sqrt(8.0)
ys[inde+3]=ycc-b/math.sqrt(8.0)
zs[inde+3]=zcc+b/math.sqrt(8.0)
types[inde+3]=3
molid0[inde+3]=8

xs[inde+4]=xcc+b/math.sqrt(8.0)
ys[inde+4]=ycc+b/math.sqrt(8.0)
zs[inde+4]=zcc-b/math.sqrt(8.0)
types[inde+4]=3
molid0[inde+4]=8

maxx=3
maxy=3
maxz=3
nt=0
xcenter=0.0
ycenter=0.0
zcenter=0.0
cellcount=0
molpercell=8
bondpercellAE=4*10
bondpercellB=4*10
bondpercell=bondpercellAE+bondpercellB

totalatoms= atompercell * maxx * maxy * maxz
xpos = [0.00] * totalatoms
ypos = [0.00] * totalatoms
zpos = [0.00] * totalatoms
typepos = [0] * totalatoms
molidcore= [0] * totalatoms

tp1count=0
tp3count=0
tp4count=0
idoftp1=[]
idoftp3=[]
idoftp4=[]
for nx in range(maxx):
  for ny in range(maxy):
    for nz in range(maxz):
      cellcount += 1
      for i in range(atompercell):
        xpos[nt]=xs[i]+nx*a
        ypos[nt]=ys[i]+ny*a
        zpos[nt]=zs[i]+nz*a
        typepos[nt]=types[i]
        if typepos[nt]==1:
          tp1count += 1
          idoftp1.append(nt)
        if typepos[nt]==3:
          tp3count += 1
          idoftp3.append(nt)
        if typepos[nt]==4:
          tp4count += 1
          idoftp4.append(nt)
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
xhi=xlo+(maxx)*a
zoff=0.0
boxx=xhi-xlo
boxy=boxx
boxz=boxx

bondnum=bondpercell*cellcount

atomtypes=4
bondtypescount=4
print("all bonds= ", bondnum)
flammps = open("TETpTET_lammpsinput_3by3.txt", "w")
flammps.write("LAMMPS input data for multiple clusters\r\n\r\n")
flammps.write("%d atoms\r\n" % nt)
flammps.write("%d bonds\r\n" % bondnum)

flammps.write("%d atom types\r\n" % atomtypes)
flammps.write("%d bond types\r\n" % bondtypescount)

flammps.write("%f %f xlo xhi\r\n" % (xlo, xhi))
flammps.write("%f %f ylo yhi\r\n" % (xlo, xhi))
flammps.write("%f %f zlo zhi\r\n" % (xlo, xhi))

flammps.write("\nMasses\n\n")
flammps.write("%d %fe-14\r\n" % (1, massE*1e14))
flammps.write("%d %fe-14\r\n" % (2, massA*1e14))
flammps.write("%d %fe-14\r\n" % (3, massB*1e14))
flammps.write("%d %fe-14\r\n" % (4, massD*1e14))


#Atoms
flammps.write("\nAtoms\r\n\n")
for i0 in range(nt):
  flammps.write("%d %d %d %f %f %f\r\n" % (i0+1, molidcore[i0], typepos[i0], xpos[i0], ypos[i0], zpos[i0]))
  
#Bonds
flammps.write("\nBonds\r\n\n")
bn=0

for i in range(cellcount):
  bn = i* bondpercell
  molstart=i*atompercell
  for j in range(4): #for tet type 1
    bondtype=1
    at_n1=molstart+j*5+1
    at_n2=molstart+j*5+2
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+1, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*5+1
    at_n2=molstart+j*5+3
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+2, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*5+1
    at_n2=molstart+j*5+4
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+3, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*5+2
    at_n2=molstart+j*5+3
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+4, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*5+2
    at_n2=molstart+j*5+4
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+5, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*5+3
    at_n2=molstart+j*5+4
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+6, bondtype,at_n1,at_n2))
    
    bondtype=2
    at_n1=molstart+j*5+1
    at_n2=molstart+j*5+5
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+7, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*5+2
    at_n2=molstart+j*5+5
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+8, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*5+3
    at_n2=molstart+j*5+5
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+9, bondtype,at_n1,at_n2))
    
    at_n1=molstart+j*5+4
    at_n2=molstart+j*5+5
    flammps.write("%d %d %d %d\r\n" % (bn+j*10+10, bondtype,at_n1,at_n2))
    
  for j in range(4): #for tet type 2
    bondtype=3
    at_n1=molstart+20+j*5+1
    at_n2=molstart+20+j*5+2
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+1, bondtype,at_n1,at_n2))
    
    at_n1=molstart+20+j*5+1
    at_n2=molstart+20+j*5+3
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+2, bondtype,at_n1,at_n2))
    
    at_n1=molstart+20+j*5+1
    at_n2=molstart+20+j*5+4
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+3, bondtype,at_n1,at_n2))
    
    at_n1=molstart+20+j*5+2
    at_n2=molstart+20+j*5+3
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+4, bondtype,at_n1,at_n2))
    
    at_n1=molstart+20+j*5+2
    at_n2=molstart+20+j*5+4
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+5, bondtype,at_n1,at_n2))
    
    at_n1=molstart+20+j*5+3
    at_n2=molstart+20+j*5+4
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+6, bondtype,at_n1,at_n2))
    
    bondtype=4
    at_n1=molstart+20+j*5+1
    at_n2=molstart+20+j*5+5
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+7, bondtype,at_n1,at_n2))
    
    at_n1=molstart+20+j*5+2
    at_n2=molstart+20+j*5+5
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+8, bondtype,at_n1,at_n2))
    
    at_n1=molstart+20+j*5+3
    at_n2=molstart+20+j*5+5
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+9, bondtype,at_n1,at_n2))
    
    at_n1=molstart+20+j*5+4
    at_n2=molstart+20+j*5+5
    flammps.write("%d %d %d %d\r\n" % (bn+bondpercellAE+j*10+10, bondtype,at_n1,at_n2))
    
flammps.close()
print("done with lammps input file")