clc;clear all;
%Schmidt Scheme for parabolic type partial differential equation
%pdefile(a,b)=y0 y0
x=1;t=0.05;
h=0.1;
k=5/1000;
c=[];
a=zeros((t/k)+1,(x/h)+1);
for j=2:1:5
  a(1,j)=2*(j-1)*h;
endfor
for j=6:1:10
  a(1,j)=2*(1-(j-1)*h);
endfor
for i=2:1:(t/k)+1
  for j=2:1:x/h
    a(i,j)=a(i-1,j)+(k/h^2)*(a(i-1,j-1)-2*a(i-1,j)+a(i-1,j+1));
  endfor
endfor
for f=0:k:t
  b=[];
  for d=0:h:1
    u=pdefile(d,f);
    b=[b u];
  endfor
  c=[c; b];
endfor
q=[0:h:x];
plot(q,a(11,:))
hold all
plot(q,c(11,:))
xlabel('spacial distance')
ylabel('time');
title("Analytical solution vs Schmidt method");
legend('numerical solution','analytical solution')

