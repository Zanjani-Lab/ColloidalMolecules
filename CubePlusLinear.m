% This code defines a superstructure of Cubic and Linear Colloidal Clusters
% Details about the structure are presented in the following reference:
% Aryana, K., Stahley, J. B., Parvez, N., Kim, K. and Zanjani, M. B. (2019) Superstructures of multielement colloidal molecules: efficient pathways to construct reconfigurable photonic and phononic crystals. Adv. Theory Simul. p. 1800198.
% Author: Zanjani research group: http://cec.miamioh.edu/zanjani-research 
% For details about the code setup, visualization in Ovito, and performing MD simulations refer to the code provided for Cubic + Octahedral units
% particle sizes are chosen for one specific scenario, they can be freely changed according to geometric constraints discussed in the article
clear;
clc;

dB=400e-7;
dBeq=dB+11e-7;
poldensity=1.05; 
massB=poldensity*pi*dB^3/6;

b=dBeq;
rEB=1.0;
dE=rEB*dB;
a=dBeq*0.5*(1+rEB+sqrt((1+rEB)^2-2.0));


dA=1.5*dE;
e=dE+11.0e-7;
dD=dB;

massE=poldensity*pi*dE^3/6;



%define (1) BD cube
inde=9;
xcc=0.0;ycc=0.0;zcc=0.0;
xs(inde+9)=xcc;ys(inde+9)=ycc;zs(inde+9)=zcc;types(inde+9)=4;molid0(inde+9)=1;
xs(inde+1)=xcc-b/2.0;ys(inde+1)=ycc-b/2.0;zs(inde+1)=zcc-b/2.0;types(inde+1)=3;molid0(inde+1)=1;
xs(inde+2)=xcc+b/2.0;ys(inde+2)=ycc-b/2.0;zs(inde+2)=zcc-b/2.0;types(inde+2)=3;molid0(inde+2)=1;
xs(inde+3)=xcc+b/2.0;ys(inde+3)=ycc+b/2.0;zs(inde+3)=zcc-b/2.0;types(inde+3)=3;molid0(inde+3)=1;
xs(inde+4)=xcc-b/2.0;ys(inde+4)=ycc+b/2.0;zs(inde+4)=zcc-b/2.0;types(inde+4)=3;molid0(inde+4)=1;

xs(inde+5)=xcc-b/2.0;ys(inde+5)=ycc-b/2.0;zs(inde+5)=zcc+b/2.0;types(inde+5)=3;molid0(inde+5)=1;
xs(inde+6)=xcc+b/2.0;ys(inde+6)=ycc-b/2.0;zs(inde+6)=zcc+b/2.0;types(inde+6)=3;molid0(inde+6)=1;
xs(inde+7)=xcc+b/2.0;ys(inde+7)=ycc+b/2.0;zs(inde+7)=zcc+b/2.0;types(inde+7)=3;molid0(inde+7)=1;
xs(inde+8)=xcc-b/2.0;ys(inde+8)=ycc+b/2.0;zs(inde+8)=zcc+b/2.0;types(inde+8)=3;molid0(inde+8)=1;


%define (3) AE linkers (linear units)
inde=0;
xcc=a;ycc=0.0;zcc=0.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=2;%A
xs(inde+1)=xcc-e/2.0;ys(inde+1)=ycc;zs(inde+1)=zcc;types(inde+1)=1;molid0(inde+1)=2;%E
xs(inde+2)=xcc+e/2.0;ys(inde+2)=ycc;zs(inde+2)=zcc;types(inde+2)=1;molid0(inde+2)=2;

inde=3;
xcc=0.0;ycc=a;zcc=0.0;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=3;%A
xs(inde+1)=xcc;ys(inde+1)=ycc-e/2.0;zs(inde+1)=zcc;types(inde+1)=1;molid0(inde+1)=3;%E
xs(inde+2)=xcc;ys(inde+2)=ycc+e/2.0;zs(inde+2)=zcc;types(inde+2)=1;molid0(inde+2)=3;

inde=6;
xcc=0.0;ycc=0.0;zcc=a;
xs(inde+3)=xcc;ys(inde+3)=ycc;zs(inde+3)=zcc;types(inde+3)=2;molid0(inde+3)=4;%A
xs(inde+1)=xcc;ys(inde+1)=ycc;zs(inde+1)=zcc-e/2.0;types(inde+1)=1;molid0(inde+1)=4;%E
xs(inde+2)=xcc;ys(inde+2)=ycc;zs(inde+2)=zcc+e/2.0;types(inde+2)=1;molid0(inde+2)=4;


maxx=3;
maxy=3;
maxz=3;
nt=0;
xcenter=0.0;ycenter=0.0;zcenter=0.0;
cellcount=0;
molpercell=4;
bondpercellAE=3*3;%linear
bondpercellB=1*20; %cube
bondpercell=bondpercellAE+bondpercellB;
atompercell=18;
for nx=0:maxx
    for ny=0:maxy
        for nz=0:maxz
            cellcount=cellcount+1;
            for i=1:atompercell
            nt=nt+1;
            xpos(nt)=xs(i)+nx*a*2;
            ypos(nt)=ys(i)+ny*a*2;
            zpos(nt)=zs(i)+nz*a*2;
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

xlo=corexmin;
xhi=xlo+(maxx+1)*2*a;
zoff=0.0;
bondnum=bondpercell*cellcount;

fid0=fopen('Seed_CubeLinear.xyz','w');
fid=fopen('LammpsInput_SeedCubeLinear.txt','w');

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
    for j=1:3
    
       
    bondtype=1; %between E particles
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+(j-1)*3+1,bondtype,molstart+(j-1)*3+1,molstart+(j-1)*3+2);
        
    bondtype=2; %between A and E particle
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+(j-1)*3+2,bondtype,molstart+(j-1)*3+1,molstart+(j-1)*3+3);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+(j-1)*3+3,bondtype,molstart+(j-1)*3+2,molstart+(j-1)*3+3);
    
    end
    
    for j=1:1
          bondtype=3; %between B particles
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+1,bondtype,molstart+9+1,molstart+9+2);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+2,bondtype,molstart+9+2,molstart+9+3);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+3,bondtype,molstart+9+3,molstart+9+4);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+4,bondtype,molstart+9+1,molstart+9+4);
    
   
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+5,bondtype,molstart+9+5,molstart+9+6);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+6,bondtype,molstart+9+6,molstart+9+7);    
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+7,bondtype,molstart+9+7,molstart+9+8);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+8,bondtype,molstart+9+5,molstart+9+8);
    
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+9,bondtype,molstart+9+1,molstart+9+5);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+10,bondtype,molstart+9+2,molstart+9+6);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+11,bondtype,molstart+9+3,molstart+9+7);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+12,bondtype,molstart+9+4,molstart+9+8);
    
    
     bondtype=4; %between B and D particle
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+13,bondtype,molstart+9+1,molstart+9+9);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+14,bondtype,molstart+9+2,molstart+9+9);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+15,bondtype,molstart+9+3,molstart+9+9);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+16,bondtype,molstart+9+4,molstart+9+9);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+17,bondtype,molstart+9+5,molstart+9+9);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+18,bondtype,molstart+9+6,molstart+9+9);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+19,bondtype,molstart+9+7,molstart+9+9);
    fprintf(fid,'%d    %d    %d    %d\r\n',bn+bondpercellAE+(j-1)*20+20,bondtype,molstart+9+8,molstart+9+9);
    
    
    end
          
end   
fclose(fid0);
fclose(fid);
disp('finished')