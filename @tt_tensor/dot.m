function [p] = dot(tt1,tt2)
%[PR]=DOT(TT1,TT2)
%Dot  product of two TT tensors
%
%   See also TT_TENSOR.
%

d=tt1.d; 
r1=tt1.r; r2=tt2.r; ps1=tt1.ps; ps2=tt2.ps;
n=tt1.n;
core1=tt1.core; core2=tt2.core;
p=1; %Vector p 
for i=1:d
  cr1=core1(ps1(i):ps1(i+1)-1); cr2=core2(ps2(i):ps2(i+1)-1);
  cr1=reshape(cr1,[r1(i),n(i)*r1(i+1)]);
p=p'*cr1; %p is now r2(i)*n(i)*r1(i+1), sum over r2(i)*n(i);
  p=reshape(p,[r2(i)*n(i),r1(i+1)]);
  cr2=reshape(cr2,[r2(i)*n(i),r2(i+1)]);
  p=p'*cr2;
end
return
end