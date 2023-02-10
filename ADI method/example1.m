clc;clear all;
%ADI method for 2D+time equations
x=0.3;y=0.3;t=0.01;
dt=0.01;dx=0.1;dy=0.1;
fx=dt/(dx^2);fy=dt/(dy^2);
t1=floor(t/dt);
x1=floor(x/dx)+1;y1=floor(y/dy)+1;
a=zeros(x1+1,y1+1,2*t1+1);
a(:,1,:)=80;a(4,:,:)=80;
a(:,:,1)=0;
for fi=2:1:2*t1+1
  if(rem(fi,2)==0)
    for o=2:1:y1
      mat=zeros(x1-1);
      mat(1,1)=-2*(fx+1);
      mat(1,2)=fx;
      g=-fy*a(o-1,2,fi-1)+2*(fy-1)*a(o,2,fi-1)-fy*a(o+1,2,fi-1)-fx*a(o,1,fi);
      b=[g];c=[];
      for i=2:1:x1-2
        mat(i,i-1)= fx;
        mat(i,i)=-2*(fx+1);
        mat(i,i+1)=fx;
        g=-fy*a(o-1,i+1,fi-1)+2*(fy-1)*a(o,i+1,fi-1)-fy*a(o+1,i+1,fi-1);
        b=[b;g];
      endfor
      mat(2,1)=fx;
      mat(2,2)=-2*(fx+1);
      g=-fy*a(o-1,3,fi-1)+(2*fy-2)*a(o,3,fi-1)-fy*a(o+1,3,fi-1)-fx*a(o,3+1,fi);
      b=[b;g];
      fs=gauss3(mat,b);
      fs
      for l=2:1:3
       a(o,l,fi)=fs(1,l-1);
      endfor
    endfor
  elseif(rem(fi,2)!=0)
    for o=2:1:x1
      mat=zeros(2);
      mat(1,1)=-2*(fy+1);
      mat(1,2)=fy;
      g=-fx*a(2,o-1,fi-1)+2*(fx-1)*a(2,o,fi-1)-fx*a(2,o+1,fi-1)-fy*a(1,o,fi);
      b=[g];c=[];
      for i=2:1:y1-2
        mat(i,i-1)= fy;
        mat(i,i)=-2*(fy+1);
        mat(i,i+1)=fy;
        g=-fx*a(i+1,o-1,fi-1)+2*(fx-1)*a(i+1,o,fi-1)-fx*a(i+1,o+1,fi-1);
        b=[b;g];
      endfor
      mat(2,1)=fy;
      mat(2,2)=-2*(fy+1);
      g=-fx*a(y1,o-1,fi-1)+(2*fx-2)*a(y1,o,fi-1)-fx*a(y1,o+1,fi-1)-fy*a(y1+1,o,fi);
      b=[b;g];
      fs=gauss3(mat,b);
      for l=2:1:y1
       a(l,o,fi)=fs(1,l-1);
      endfor
    endfor
  endif
endfor
