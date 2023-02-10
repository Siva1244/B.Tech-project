function [f]=pdefile(a,b)
k=0;
for j=1:1:10  
  k=k+sum((1/(j.^2)).*sin(j*pi/2).*sin(j.*pi.*a).*exp(-j.^2*pi^2.*b));
endfor
f=(8/pi^2)*k;
endfunction  