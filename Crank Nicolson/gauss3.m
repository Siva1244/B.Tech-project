%matrix reduction(row reduced echolon form)
function f=gauss3(a,mat)
% Siva on 5th Jan,2018 
b=size(a,1);
c=size(a,2);
a=[a mat];
f=zeros(1,b);
for j=1:1:b
  for k=j:1:b-1
    if(a(j,j)==0) % check for pivot element
     f(1,:)=a(j,:);
     a(j,:)=a(k+1,:);
     a(k+1,:)=f(1,:);
    endif
  endfor
  if(a(j,j)!=0)
    a(j,:)=a(j,:)/a(j,j);
   for i=j+1:1:b
        a(i,:)=a(i,:)-a(i,j)*a(j,:);
   endfor
  endif
endfor
f(1,b)=a(b,c+1);
for i=b-1:-1:1
  s=a(i,c+1);
  r=b;
  for e=c:-1:i+1
    s=s-a(i,e)*f(1,r);
    r=r-1;
  endfor
  f(1,r)=s;
endfor
endfunction