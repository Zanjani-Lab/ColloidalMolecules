% This MATLAB code defines a superstructure of Tetrahedral and Linear Colloidal Clusters
% Details about the structure are presented in the following reference:
% Aryana, K., Stahley, J. B., Parvez, N., Kim, K. and Zanjani, M. B. (2019) Superstructures of multielement colloidal molecules: efficient pathways to construct reconfigurable photonic and phononic crystals. Adv. Theory Simul. p. 1800198.
% Author: Zanjani research group: http://cec.miamioh.edu/zanjani-research 
% For details about the code setup, visualization in Ovito, and performing MD simulations refer to the code provided for Cubic + Octahedral units
% particle sizes are chosen for one specific scenario, they can be freely changed according to geometric constraints discussed in the article
clear;
clc;

dB=400e-7;
dBeq=dB+11e-7;
g=11.0e-07;
poldensity=1.05; 
massB=poldensity*pi*dB^3/6;
dD=0.66*dB;
dDeq=dD+11.0e-7;

s=(g^2+dB*dD+g*dB+g*dD)/((sqrt(3)+1/sqrt(3))*dB-(sqrt(3)-1/sqrt(3))*dD+2*g/sqrt(3));
dE=2*sqrt(3)*s-dD-2*g;
rEB=dE/dB;
e=(dE+11.0e-7);
a=8*s+4*e/sqrt(3)
dA=1.5*dE;

b=sqrt(8.0/3.0)*(dDeq+dBeq)/2; 
%define BD cubes
inde0=48;
inde=inde0;
xcc=0.0;ycc=0.0;zcc=0.0;
xs(inde+5)=xcc;ys(inde+5)=ycc;zs(inde+5)=zcc;types(inde+5)=4;molid0(inde+5)=1; %B
xs(inde+1)=xcc+b/sqrt(8.0);ys(inde+1)=ycc-b/sqrt(8.0);zs(inde+1)=zcc+b/sqrt(8.0);types(inde+1)=3;molid0(inde+1)=1;%D
xs(inde+2)=xcc-b/sqrt(8.0);ys(inde+2)=ycc+b/sqrt(8.0);zs(inde+2)=zcc+b/sqrt(8.0);types(inde+2)=3;molid0(inde+2)=1;
xs(inde+3)=xcc+b/sqrt(8.0);ys(inde+3)=ycc+b/sqrt(8.0);zs(inde+3)=zcc-b/sqrt(8.0);types(inde+3)=3;molid0(inde+3)=1;
xs(inde+4)=xcc-b/sqrt(8.0);ys(inde+4)=ycc-b/sqrt(8.0);zs(inde+4)=zcc-b/sqrt(8.0);types(inde+4)=3;molid0(inde+4)=1;

inde=inde0+5;
xcc=0.5*a;ycc=0.5*a;zcc=0.0;
xs(inde+5)=xcc;ys(inde+5)=ycc;zs(inde+5)=zcc;types(inde+5)=4;molid0(inde+5)=2;
xs(inde+1)=xcc+b/sqrt(8.0);ys(inde+1)=ycc-b/sqrt(8.0);zs(inde+1)=zcc+b/sqrt(8.0);types(inde+1)=3;molid0(inde+1)=2;
xs(inde+2)=xcc-b/sqrt(8.0);ys(inde+2)=ycc+b/sqrt(8.0);zs(inde+2)=zcc+b/sqrt(8.0);types(inde+2)=3;molid0(inde+2)=2;
xs(inde+3)=xcc+b/sqrt(8.0);ys(inde+3)=ycc+b/sqrt(8.0);zs(inde+3)=zcc-b/sqrt(8.0);types(inde+3)=3;molid0(inde+3)=2;
xs(inde+4)=xcc-b/sqrt(8.0);ys(inde+4)=ycc-b/sqrt(8.0);zs(inde+4)=zcc-b/sqrt(8.0);types(inde+4)=3;molid0(inde+4)=2;

inde=inde0+10;
xcc=0.0;ycc=0.5*a;zcc=0.5*a;
xs(inde+5)=xcc;ys(inde+5)=ycc;zs(inde+5)=zcc;types(inde+5)=4;molid0(inde+5)=3;
xs(inde+1)=xcc+b/sqrt(8.0);ys(inde+1)=ycc-b/sqrt(8.0);zs(inde+1)=zcc+b/sqrt(8.0);types(inde+1)=3;molid0(inde+1)=3;
xs(inde+2)=xcc-b/sqrt(8.0);ys(inde+2)=ycc+b/sqrt(8.0);zs(inde+2)=zcc+b/sqrt(8.0);types(inde+2)=3;molid0(inde+2)=3;
xs(inde+3)=xcc+b/sqrt(8.0);ys(inde+3)=ycc+b/sqrt(8.0);zs(inde+3)=zcc-b/sqrt(8.0);types(inde+3)=3;molid0(inde+3)=3;
xs(inde+4)=xcc-b/sqrt(8.0);ys(inde+4)=ycc-b/sqrt(8.0);zs(inde+4)=zcc-b/sqrt(8.0);types(inde+4)=3;molid0(inde+4)=3;

inde=inde0+15;
xcc=0.5*a;ycc=0.0;zcc=0.5*a;
xs(inde+5)=xcc;ys(inde+5)=ycc;zs(inde+5)=zcc;types(inde+5)=4;molid0(inde+5)=4;
xs(inde+1)=xcc+b/sqrt(8.0);ys(inde+1)=ycc-b/sqrt(8.0);zs(inde+1)=zcc+b/sqrt(8.0);types(inde+1)=3;molid0(inde+1)=4;
xs(inde+2)=xcc-b/sqrt(8.0);ys(inde+2)=ycc+b/sqrt(8.0);zs(inde+2)=zcc+b/sqrt(8.0);types(inde+2)=3;molid0(inde+2)=4;
xs(inde+3)=xcc+b/sqrt(8.0);ys(inde+3)=ycc+b/sqrt(8.0);zs(inde+3)=zcc-b/sqrt(8.0);types(inde+3)=3;molid0(inde+3)=4;
xs(inde+4)=xcc-b/sqrt(8.0);ys(inde+4)=ycc-b/sqrt(8.0);zs(inde+4)=zcc-b/sqrt(8.0);types(inde+4)=3;molid0(inde+4)=4;

inde=inde0+20;
xcc=0.25*a;ycc=0.25*a;zcc=0.25*a;
xs(inde+5)=xcc;ys(inde+5)=ycc;zs(inde+5)=zcc;types(inde+5)=4;molid0(inde+5)=5;
xs(inde+1)=xcc+b/sqrt(8.0);ys(inde+1)=ycc+b/sqrt(8.0);zs(inde+1)=zcc+b/sqrt(8.0);types(inde+1)=3;molid0(inde+1)=5;
xs(inde+2)=xcc+b/sqrt(8.0);ys(inde+2)=ycc-b/sqrt(8.0);zs(inde+2)=zcc-b/sqrt(8.0);types(inde+2)=3;molid0(inde+2)=5;
xs(inde+3)=xcc-b/sqrt(8.0);ys(inde+3)=ycc+b/sqrt(8.0);zs(inde+3)=zcc-b/sqrt(8.0);types(inde+3)=3;molid0(inde+3)=5;
xs(inde+4)=xcc-b/sqrt(8.0);ys(inde+4)=ycc-b/sqrt(8.0);zs(inde+4)=zcc+b/sqrt(8.0);types(inde+4)=3;molid0(inde+4)=5;

inde=inde0+25;
xcc=0.75*a;ycc=0.75*a;zcc=0.25*a;
xs(inde+5)=xcc;ys(inde+5)=ycc;zs(inde+5)=zcc;types(inde+5)=4;molid0(inde+5)=6;
xs(inde+1)=xcc+b/sqrt(8.0);ys(inde+1)=ycc+b/sqrt(8.0);zs(inde+1)=zcc+b/sqrt(8.0);types(inde+1)=3;molid0(inde+1)=6;
xs(inde+2)=xcc+b/sqrt(8.0);ys(inde+2)=ycc-b/sqrt(8.0);zs(inde+2)=zcc-b/sqrt(8.0);types(inde+2)=3;molid0(inde+2)=6;
xs(inde+3)=xcc-b/sqrt(8.0);ys(inde+3)=ycc+b/sqrt(8.0);zs(inde+3)=zcc-b/sqrt(8.0);types(inde+3)=3;molid0(inde+3)=6;
xs(inde+4)=xcc-b/sqrt(8.0);ys(inde+4)=ycc-b/sqrt(8.0);zs(inde+4)=zcc+b/sqrt(8.0);types(inde+4)=3;molid0(inde+4)=6;

inde=inde0+30;
xcc=0.25*a;ycc=0.75*a;zcc=0.75*a;
xs(inde+5)=xcc;ys(inde+5)=ycc;zs(inde+5)=zcc;types(inde+5)=4;molid0(inde+5)=7;
xs(inde+1)=xcc+b/sqrt(8.0);ys(inde+1)=ycc+b/sqrt(8.0);zs(inde+1)=zcc+b/sqrt(8.0);types(inde+1)=3;molid0(inde+1)=7;
xs(inde+2)=xcc+b/sqrt(8.0);ys(inde+2)=ycc-b/sqrt(8.0);zs(inde+2)=zcc-b/sqrt(8.0);types(inde+2)=3;molid0(inde+2)=7;
xs(inde+3)=xcc-b/sqrt(8.0);ys(inde+3)=ycc+b/sqrt(8.0);zs(inde+3)=zcc-b/sqrt(8.0);types(inde+3)=3;molid0(inde+3)=7;
xs(inde+4)=xcc-b/sqrt(8.0);ys(inde+4)=ycc-b/sqrt(8.0);zs(inde+4)=zcc+b/sqrt(8.0);types(inde+4)=3;molid0(inde+4)=7;

inde=inde0+35;
xcc=0.75*a;ycc=0.25*a;zcc=0.75*a;
xs(inde+5)=xcc;ys(inde+5)=ycc;zs(inde+5)=zcc;types(inde+5)=4;molid0(inde+5)=8;
xs(inde+1)=xcc+b/sqrt(8.0);ys(inde+1)=ycc+b/sqrt(8.0);zs(inde+1)=zcc+b/sqrt(8.0);types(inde+1)=3;molid0(inde+1)=8;
xs(inde+2)=xcc+b/sqrt(8.0);ys(inde+2)=ycc-b/sqrt(8.0);zs(inde+2)=zcc-b/sqrt(8.0);types(inde+2)=3;molid0(inde+2)=8;
xs(inde+3)=xcc-b/sqrt(8.0);ys(inde+3)=ycc+b/sqrt(8.0);zs(inde+3)=zcc-b/sqrt(8.0);types(inde+3)=3;molid0(inde+3)=8;
xs(inde+4)=xcc-b/sqrt(8.0);ys(inde+4)=ycc-b/sqrt(8.0);zs(inde+4)=zcc+b/sqrt(8.0);types(inde+4)=3;molid0(inde+4)=8;

% define AE linear

inde=0;
xcc=a/4.0-a/8.0;ycc=a/4.0-a/8.0;zcc=a/4.0-a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=9;%A
xs(inde+1)=xcc-0.5*e/sqrt(3.0);ys(inde+1)=ycc-0.5*e/sqrt(3.0);zs(inde+1)=zcc-0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=9;%E
xs(inde+2)=xcc+0.5*e/sqrt(3.0);ys(inde+2)=ycc+0.5*e/sqrt(3.0);zs(inde+2)=zcc+0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=9;

inde=3;
xcc=a/4.0+a/8.0;ycc=a/4.0-a/8.0;zcc=a/4.0+a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=10;%A
xs(inde+1)=xcc+0.5*e/sqrt(3.0);ys(inde+1)=ycc-0.5*e/sqrt(3.0);zs(inde+1)=zcc+0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=10;%E
xs(inde+2)=xcc-0.5*e/sqrt(3.0);ys(inde+2)=ycc+0.5*e/sqrt(3.0);zs(inde+2)=zcc-0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=10;

inde=6;
xcc=a/4.0-a/8.0;ycc=a/4.0+a/8.0;zcc=a/4.0+a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=11;%A
xs(inde+1)=xcc-0.5*e/sqrt(3.0);ys(inde+1)=ycc+0.5*e/sqrt(3.0);zs(inde+1)=zcc+0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=11;%E
xs(inde+2)=xcc+0.5*e/sqrt(3.0);ys(inde+2)=ycc-0.5*e/sqrt(3.0);zs(inde+2)=zcc-0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=11;

inde=9;
xcc=a/4.0+a/8.0;ycc=a/4.0+a/8.0;zcc=a/4.0-a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=12;%A
xs(inde+1)=xcc+0.5*e/sqrt(3.0);ys(inde+1)=ycc+0.5*e/sqrt(3.0);zs(inde+1)=zcc-0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=12;%E
xs(inde+2)=xcc-0.5*e/sqrt(3.0);ys(inde+2)=ycc-0.5*e/sqrt(3.0);zs(inde+2)=zcc+0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=12;

inde=12;
xcc=3*a/4.0-a/8.0;ycc=3*a/4.0-a/8.0;zcc=a/4.0-a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=13;%A
xs(inde+1)=xcc-0.5*e/sqrt(3.0);ys(inde+1)=ycc-0.5*e/sqrt(3.0);zs(inde+1)=zcc-0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=13;%E
xs(inde+2)=xcc+0.5*e/sqrt(3.0);ys(inde+2)=ycc+0.5*e/sqrt(3.0);zs(inde+2)=zcc+0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=13;

inde=15;
xcc=3*a/4.0+a/8.0;ycc=3*a/4.0-a/8.0;zcc=a/4.0+a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=14;%A
xs(inde+1)=xcc+0.5*e/sqrt(3.0);ys(inde+1)=ycc-0.5*e/sqrt(3.0);zs(inde+1)=zcc+0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=14;%E
xs(inde+2)=xcc-0.5*e/sqrt(3.0);ys(inde+2)=ycc+0.5*e/sqrt(3.0);zs(inde+2)=zcc-0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=14;

inde=18;
xcc=3*a/4.0-a/8.0;ycc=3*a/4.0+a/8.0;zcc=a/4.0+a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=15;%A
xs(inde+1)=xcc-0.5*e/sqrt(3.0);ys(inde+1)=ycc+0.5*e/sqrt(3.0);zs(inde+1)=zcc+0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=15;%E
xs(inde+2)=xcc+0.5*e/sqrt(3.0);ys(inde+2)=ycc-0.5*e/sqrt(3.0);zs(inde+2)=zcc-0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=15;

inde=21;
xcc=3*a/4.0+a/8.0;ycc=3*a/4.0+a/8.0;zcc=a/4.0-a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=16;%A
xs(inde+1)=xcc+0.5*e/sqrt(3.0);ys(inde+1)=ycc+0.5*e/sqrt(3.0);zs(inde+1)=zcc-0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=16;%E
xs(inde+2)=xcc-0.5*e/sqrt(3.0);ys(inde+2)=ycc-0.5*e/sqrt(3.0);zs(inde+2)=zcc+0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=16;

inde=24;
xcc=a/4.0-a/8.0;ycc=3*a/4.0-a/8.0;zcc=3*a/4.0-a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=17;%A
xs(inde+1)=xcc-0.5*e/sqrt(3.0);ys(inde+1)=ycc-0.5*e/sqrt(3.0);zs(inde+1)=zcc-0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=17;%E
xs(inde+2)=xcc+0.5*e/sqrt(3.0);ys(inde+2)=ycc+0.5*e/sqrt(3.0);zs(inde+2)=zcc+0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=17;

inde=27;
xcc=a/4.0+a/8.0;ycc=3*a/4.0-a/8.0;zcc=3*a/4.0+a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=18;%A
xs(inde+1)=xcc+0.5*e/sqrt(3.0);ys(inde+1)=ycc-0.5*e/sqrt(3.0);zs(inde+1)=zcc+0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=18;%E
xs(inde+2)=xcc-0.5*e/sqrt(3.0);ys(inde+2)=ycc+0.5*e/sqrt(3.0);zs(inde+2)=zcc-0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=18;

inde=30;
xcc=a/4.0-a/8.0;ycc=3*a/4.0+a/8.0;zcc=3*a/4.0+a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=19;%A
xs(inde+1)=xcc-0.5*e/sqrt(3.0);ys(inde+1)=ycc+0.5*e/sqrt(3.0);zs(inde+1)=zcc+0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=19;%E
xs(inde+2)=xcc+0.5*e/sqrt(3.0);ys(inde+2)=ycc-0.5*e/sqrt(3.0);zs(inde+2)=zcc-0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=19;

inde=33;
xcc=a/4.0+a/8.0;ycc=3*a/4.0+a/8.0;zcc=3*a/4.0-a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=20;%A
xs(inde+1)=xcc+0.5*e/sqrt(3.0);ys(inde+1)=ycc+0.5*e/sqrt(3.0);zs(inde+1)=zcc-0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=20;%E
xs(inde+2)=xcc-0.5*e/sqrt(3.0);ys(inde+2)=ycc-0.5*e/sqrt(3.0);zs(inde+2)=zcc+0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=20;

inde=36;
xcc=3*a/4.0-a/8.0;ycc=a/4.0-a/8.0;zcc=3*a/4.0-a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=21;%A
xs(inde+1)=xcc-0.5*e/sqrt(3.0);ys(inde+1)=ycc-0.5*e/sqrt(3.0);zs(inde+1)=zcc-0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=21;%E
xs(inde+2)=xcc+0.5*e/sqrt(3.0);ys(inde+2)=ycc+0.5*e/sqrt(3.0);zs(inde+2)=zcc+0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=21;

inde=39;
xcc=3*a/4.0+a/8.0;ycc=a/4.0-a/8.0;zcc=3*a/4.0+a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=22;%A
xs(inde+1)=xcc+0.5*e/sqrt(3.0);ys(inde+1)=ycc-0.5*e/sqrt(3.0);zs(inde+1)=zcc+0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=22;%E
xs(inde+2)=xcc-0.5*e/sqrt(3.0);ys(inde+2)=ycc+0.5*e/sqrt(3.0);zs(inde+2)=zcc-0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=22;

inde=42;
xcc=3*a/4.0-a/8.0;ycc=a/4.0+a/8.0;zcc=3*a/4.0+a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=23;%A
xs(inde+1)=xcc-0.5*e/sqrt(3.0);ys(inde+1)=ycc+0.5*e/sqrt(3.0);zs(inde+1)=zcc+0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=23;%E
xs(inde+2)=xcc+0.5*e/sqrt(3.0);ys(inde+2)=ycc-0.5*e/sqrt(3.0);zs(inde+2)=zcc-0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=23;

inde=45;
xcc=3*a/4.0+a/8.0;ycc=a/4.0+a/8.0;zcc=3*a/4.0-a/8.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=24;%A
xs(inde+1)=xcc+0.5*e/sqrt(3.0);ys(inde+1)=ycc+0.5*e/sqrt(3.0);zs(inde+1)=zcc-0.5*e/sqrt(3.0);types(inde+1)=1;molid0(inde+1)=24;%E
xs(inde+2)=xcc-0.5*e/sqrt(3.0);ys(inde+2)=ycc-0.5*e/sqrt(3.0);zs(inde+2)=zcc+0.5*e/sqrt(3.0);types(inde+2)=1;molid0(inde+2)=24;


maxx=3;
maxy=3;
maxz=3;
nt=0;
xcenter=0.0;ycenter=0.0;zcenter=0.0;
cellcount=0;
molpercell=24;
bondpercellAE=16*3;%AE
bondpercellB=8*10; %TET
bondpercell=bondpercellAE+bondpercellB;
atompercell=40+16*3;
for nx=0:maxx
    for ny=0:maxy
        for nz=0:maxz
            cellcount=cellcount+1;
            for i=1:atompercell
            nt=nt+1;
            xpos(nt)=xs(i)+nx*a;
            ypos(nt)=ys(i)+ny*a;
            zpos(nt)=zs(i)+nz*a;
            typepos(nt)=types(i);
            xcenter=xcenter+xpos(nt);
            ycenter=ycenter+ypos(nt);
            zcenter=zcenter+zpos(nt);
            molidcore(nt)=(cellcount-1)*molpercell+molid0(i);
            end
        end
    end
end
molseed=cellcount*molpercell;
xcenter=xcenter/double(nt);
ycenter=ycenter/double(nt);
zcenter=zcenter/double(nt);

for i=1:nt
   xpos(i)=xpos(i)-xcenter;
   ypos(i)=ypos(i)-ycenter;
   zpos(i)=zpos(i)-zcenter;
end

corezmin=50000;
corezmax=-50000;
coreymin=50000;
coreymax=-50000;
corexmin=50000;
corexmax=-50000;

for j=1:nt
    if(zpos(j)<corezmin)
        corezmin=zpos(j);
    end
    if(zpos(j)>corezmax)
        corezmax=zpos(j);
    end
    if(xpos(j)<corexmin)
        corexmin=xpos(j);
    end
    if(xpos(j)>corexmax)
        corexmax=xpos(j);
    end
    if(ypos(j)<coreymin)
        coreymin=ypos(j);
    end
    if(ypos(j)>coreymax)
        coreymax=ypos(j);
    end
end

fid00=fopen('SeedTETDlinkNbbwED.xyz','w');

fprintf(fid00,'%d \r\n',nt);
fprintf(fid00,'Atoms\n');

for i=1:nt
    
            fprintf(fid00,'%d      %E    %E    %E\r\n',typepos(i),1e5*xpos(i),1e5*ypos(i),1e5*zpos(i));
end

fclose(fid00);

xlo=corexmin;
xhi=xlo+(maxx+1)*a;
zoff=0.0;
bondnum=bondpercell*cellcount;

fid0=fopen('TetrahedraplusLinear.xyz','w');
fid=fopen('LammpsInput_TetrahedraPlusLinear.txt','w');

massE=poldensity*pi*dE^3/6;
massB=poldensity*pi*dB^3/6;
massA=poldensity*pi*dA^3/6;
massD=poldensity*pi*dD^3/6;

fprintf(fid0,'%d \r\n',nt);
fprintf(fid0,'Atoms\n');

fprintf(fid,'LAMMPS input data for multiple clusters\r\n\r\n');
fprintf(fid,'%d atoms\r\n',nt);
fprintf(fid,'%d bonds\r\n\n',bondnum);
fprintf(fid, '4 atom types\r\n');
fprintf(fid, '4 bond types\r\n\n');
fprintf(fid,'%f %f xlo xhi\r\n',xlo-zoff,xhi+zoff);
fprintf(fid,'%f %f ylo yhi\r\n',xlo-zoff,xhi+zoff);
fprintf(fid,'%f %f zlo zhi\r\n\n',xlo-zoff,xhi+zoff);
fprintf(fid,'Masses\n');
fprintf(fid,'\n');
fprintf(fid,'%d %E \r\n',1, massE);
fprintf(fid,'%d %E \r\n',2, massA);
fprintf(fid,'%d %E \r\n',3, massB);
fprintf(fid,'%d %E \r\n',4, massD);

fprintf(fid,'\n');

%atoms position
fprintf(fid,'Atoms\r\n\n');
k=1;
atomindx=0;
Aindx=0;
mlid=0;
for i=1:nt
    atomindx=atomindx+1;
    k=molidcore(i);
    fprintf(fid,'%d     %d     %d    %E    %E    %E\r\n',atomindx,k,typepos(i),xpos(i),ypos(i),zpos(i));
    fprintf(fid0,'%d      %E    %E    %E\r\n',typepos(i),1e5*xpos(i),1e5*ypos(i),1e5*zpos(i));
end

%bonds definition
fprintf(fid,'\nBonds\r\n\n');
bn=0;
for i=1:cellcount
    bn=(i-1)*bondpercell;
    molstart=(i-1)*atompercell;
    for j=1:16
    
       
    bondtype=1; %between E particles
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+(j-1)*3+1,bondtype,molstart+(j-1)*3+1,molstart+(j-1)*3+2);
        
    bondtype=2; %between A and E particle
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+(j-1)*3+2,bondtype,molstart+(j-1)*3+1,molstart+(j-1)*3+3);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+(j-1)*3+3,bondtype,molstart+(j-1)*3+2,molstart+(j-1)*3+3);
    
    end
    
    for j=1:8
          bondtype=3; %between B particles
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+1,bondtype,molstart+(j-1)*5+48+1,molstart+(j-1)*5+48+2);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+2,bondtype,molstart+(j-1)*5+48+1,molstart+(j-1)*5+48+3);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+3,bondtype,molstart+(j-1)*5+48+1,molstart+(j-1)*5+48+4);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+4,bondtype,molstart+(j-1)*5+48+2,molstart+(j-1)*5+48+3);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+5,bondtype,molstart+(j-1)*5+48+2,molstart+(j-1)*5+48+4);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+6,bondtype,molstart+(j-1)*5+48+3,molstart+(j-1)*5+48+4);    
    
    bondtype=4; %between B and D particles
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+7,bondtype,molstart+(j-1)*5+48+1,molstart+(j-1)*5+48+5);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+8,bondtype,molstart+(j-1)*5+48+2,molstart+(j-1)*5+48+5);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+9,bondtype,molstart+(j-1)*5+48+3,molstart+(j-1)*5+48+5);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*10+10,bondtype,molstart+(j-1)*5+48+4,molstart+(j-1)*5+48+5);
    end
          
end

fclose(fid0);
fclose(fid);
disp('finished')
