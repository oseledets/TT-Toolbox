function [p] = dot(tt1,tt2)
%[PR]=DOT(TT1,TT2)
%Dot  product of two TT tensors
%
%

[tt1,rv1]=qr(tt1, 'lr');
[tt2,rv2]=qr(tt2, 'lr');

d=tt1.d; 
r1=tt1.r; r2=tt2.r; ps1=tt1.ps; ps2=tt2.ps;
n=tt1.n;
core1=tt1.core; core2=tt2.core;
%ps is always r1(i-1)xr; but if there is a hanging thing? what to do?
%That means, I define a dot product of two "hanging" tensors as a matrix...
%brr.... And if it is hanging on the right? 
% 
% p=ones(r1(1),r2(1)); % Over r1(1) and r2(1) there is not summation blin.
% %So, if we just sum over n(1) separatedly and leave a strange QxR1(I)xR2(I)
% %matrix... 
% for i=1:d
%   cr1=core1(ps1(i):ps1(i+1)-1);
%   cr2=core2(ps2(i):ps2(i+1)-1);
%   p=reshape(p,[r1(i),r2(i)]);
%   cr2=reshape(cr2,[r2(i),numel(cr2)/r2(i)]);
%   p=p*cr2; %p is Q*r1(i)xn(i)xr2(i+1);
%   cr1=reshape(cr1,[r1(i)*n(i),numel(cr1)/(r1(i)*n(i))]);
%   p=reshape(p,[r1(i)*n(i),numel(p)/(r1(i)*n(i))]);
%   p=cr1'*p;
% end
p=1;
for i=1:d
  cr1=core1(ps1(i):ps1(i+1)-1);
  cr2=core2(ps2(i):ps2(i+1)-1);
  cr1=reshape(cr1,[r1(i),n(i)*r1(i+1)]);
  p=p'*cr1; %p is now r2(i)*n(i)*r1(i+1), sum over r2(i)*n(i);
  p=reshape(p,[r2(i)*n(i),r1(i+1)]);
  cr2=reshape(cr2,[r2(i)*n(i),r2(i+1)]);
  p=p'*cr2;
end
p = (rv1.')*p;
p = p*rv2;
return
end  
