clc;clear all;
%ADI method for 2D+time equations
x=0.3;y=0.3;t=0.03;
dt=0.01;dx=0.1;dy=0.1;
fx=dt/(dx^2);fy=dt/(dy^2);
a=zeros(y/dy+1,x/dx+1,2*(t/dt)+1)
a(:,1,:)=80;a(1,:,:)=80;
a1=x/dx-1;a2=y/dy-1;
for fi=2:1:2*(t/dt)+1
  if(rem(fi,2)==0)
    for o=2:1:y/dy
      mat=eye(2)
      mat(1,1)=-2*(fx+1);
      mat(1,2)=fx;
      g=-fy*a(o-1,2,fi-1)+2*(fy-1)*a(o,2,fi-1)-fy*a(o+1,2,fi-1)-fx*a(o,1,fi);
      b=[g];c=[];
      for i=2:1:(x/dx)-2
        mat(i,i-1)= fx;
        mat(i,i)=-2*(fx+1);
        mat(i,i+1)=fx;
        g=-fy*a(o-1,i+1,fi-1)+2*(fy-1)*a(o,i+1,fi-1)-fy*a(o+1,i+1,fi-1);
        b=[b;g];
      endfor
      mat(x/dx-1,(x/dx)-2)=fx;
      mat(x/dx-1,x/dx-1)=-2*(fx+1);
      g=-fy*a(o-1,x/dx,fi-1)+(2*fy-2)*a(o,x/dx,fi-1)-fy*a(o+1,x/dx,fi-1)-fx*a(o,x/dx+1,fi);
      b=[b;g];
      mat,b
      fs=gauss3(mat,b);
      for l=2:1:x/dx
       a(o,l,fi)=fs(1,l-1);
      endfor
    endfor
  elseif(rem(fi,2)!=0)
    for o=2:1:x/dx
      mat=zeros(y/dy-1);
      mat(1,1)=-2*(fy+1);
      mat(1,2)=fy;
      g=-fx*a(2,o-1,fi-1)+2*(fx-1)*a(2,o,fi-1)-fx*a(2,o+1,fi-1)-fy*a(1,o,fi);
      b=[g];c=[];
      for i=2:1:(y/dy)-2
        mat(i,i-1)= fy;
        mat(i,i)=-2*(fy+1);
        mat(i,i+1)=fy;
        g=-fx*a(i+1,o-1,fi-1)+2*(fx-1)*a(i+1,o,fi-1)-fx*a(i+1,o+1,fi-1);
        b=[b;g];
      endfor
      mat(y/dy-1,(y/dy)-2)=fy;
      mat(y/dy-1,y/dy-1)=-2*(fy+1);
      g=-fx*a(y/dy,o-1,fi-1)+(2*fx-2)*a(y/dy,o,fi-1)-fx*a(y/dy,o+1,fi-1)-fy*a(y/dy+1,o,fi);
      b=[b;g];
      fs=gauss3(mat,b);
      for l=2:1:y/dy
       a(l,o,fi)=fs(1,l-1);
      endfor
    endfor
  endif
endfor