function [tt] = full_to_ttr(b,ranks)
%[TT]=FULL_TO_TTR(B,RANKS)
%This subroutine gets a full d-dimensional tensor B as an input,
%accuracy RANKS are given, and TT approximation with 
%guaranteed precision ||TT-B||_EPS <= EPS*||B||_F is computed, with
%almost optimal ranks. The norm of the cores are balanced for
%enhanced stability
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
arr=b;
d=size(size(arr),2); %Dimensionality
tt=cell(d,1);
%Gradually compress modes
sz=size(arr);
nl=sz(1);
nr=1;
for i=2:d
    nr = nr*sz(i);
end
%First mode is somewhat different
nrm=norm(reshape(b,[nl*nr,1])); %For compression
%Estimate how much norm is needed for a single core
nrm1=nrm^(1.0/d);
if ( nrm < 1e-100 ) %Too small
   for i=2:d-1
     tt{i}=zeros(2,0,0);
   end
   tt{1}=zeros(2,0); tt{d}=zeros(2,0);
end
b=reshape(b,[nl,nr]);
[u0,s0,v0]=svd(b,'econ');
r=min(min(ranks(1),nl),nr);
u0=u0(:,1:r);
s0=s0(1:r,1:r);
s0=s0./nrm1;
b=v0(:,1:r)*s0;
b=b';
tt{1}=u0*nrm1;
for i=2:d-1
    nl=sz(i);
    nr=nr/sz(i);
    b=reshape(b,[r*nl,nr]);
    [u0,s0,v0]=svd(b,'econ');
    r1=min(min(ranks(i),r*nl),nr);
    u0=u0(:,1:r1).*nrm1; %u0 is rxnlxr1
    s0=s0(1:r1,1:r1)./nrm1;
    b=(v0(:,1:r1)*s0);
    b=b'; %B is r1xn2xn3
    tt{i}=reshape(u0,[r,nl,r1]);
    tt{i}=permute(tt{i},[2,1,3]);   
    r=r1;
end
tt{d}=b';
return
end
