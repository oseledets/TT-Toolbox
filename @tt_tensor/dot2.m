function [nrm]=dot2(tt1,tt2)
%Computes the scalar product of tt1 and tt2 as a product of numbers
%   [rm]=DOT2(TT1,TT2) --- rm is a vector of length d
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

d=tt1.d;
n=tt1.n;
r1=tt1.r;
r2=tt2.r;
ps1=tt1.ps;
ps2=tt2.ps;
core1=tt1.core;
core2=tt2.core;

nrm=zeros(d,1); %The norm array
p=1; %Vector p is r1(i)xr2(i)
for i=1:d
   cr1=core1(ps1(i):ps1(i+1)-1);
   cr1=reshape(cr1,[r1(i),n(i)*r1(i+1)]);
   p=p'*cr1; %p is now r2(i)*n(i)*r1(i+1), sum over r2(i)*n(i);
  p=reshape(p,[r2(i)*n(i),r1(i+1)]);
   cr2=core2(ps2(i):ps2(i+1)-1);
   cr2=reshape(cr2,[r2(i)*n(i),r2(i+1)]);
   p=p'*cr2;
   nrm(i)=norm(p,'fro');
   p=p./nrm(i);
end