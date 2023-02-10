clc;clear all;
%Crank Nicolson method for parabolic type differential equations
x=1;h=0.01;k=0.05;t=0.1;
mat1=zeros(t/k+1,x/h+1);
for i=1:1:x/h+1
  mat1(1,i)=2*sin(2*pi*(i-1)*h);
endfor
for o=2:1:3
mat=zeros(x/h);
r=k/(h^2);
mat(1,1)=(32+2*r);
mat(1,2)=-r;
g=r*mat1(o-1,1)+(32-2*r)*mat1(o-1,2)+r*mat1(o-1,3)+r*mat1(o,1);
b=[g];c=[];
for i=2:1:(x/h)-1
  mat(i,i-1)=-r;
  mat(i,i)=(32+2*r);
  mat(i,i+1)=-r;
  g=r*mat1(o-1,i-1)+(32-2*r)*mat1(o-1,i)+r*mat1(o-1,i+1);
  b=[b;g];
endfor
mat(x/h,(x/h)-1)=-r;
mat(x/h,x/h)=(32+2*r);
g=r*mat1(o-1,x/h-1)+(32-2*r)*mat1(o-1,x/h)+r*mat1(o-1,x/h+1)+r*mat1(o,x/h+1);
b=[b;g];
f=gauss3(mat,b);
for l=2:1:x/h
  mat1(o,l)=f(1,l-1);
endfor
endfor
for f1=0:k:t
  b1=[];
  for d1=0:0.01:x
    u=pdepro2(d1,f1);
    b1=[b1 u];
  endfor
  c=[c; b1];
endfor
q=[0:h:x];
plot(q,mat1((t/k)+1,:))
hold all
q=[0:0.01:x];
plot(q,c((t/k)+1,:))
mat1;