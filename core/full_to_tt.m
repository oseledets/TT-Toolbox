function [tt] = full_to_tt(b,eps,max_r)
%[TT]=FULL_TO_TT(B,EPS)
%This subroutine gets a full d-dimensional tensor B as an input,
%accuracy EPS is given, and TT approximation with 
%guaranteed precision ||TT-B||_EPS <= EPS*||B||_F is computed, with
%almost optimal ranks. The norm of the cores are balanced for
%enhanced stability
%
%
%
% TT Toolbox 1.1, 2009
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

exists_max_r = 1;
if ((nargin<3)||(isempty(max_r)))
    exists_max_r=0;
end;

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
% nrm=1;
%Estimate how much norm is needed for a single core
nrm1=nrm^(1.0/d);
% nrm1=1;
if ( nrm < 1e-307 ) %Too small
   for i=2:d-1
     tt{i}=zeros(sz(i),1,1);
   end
   tt{1}=zeros(sz(1),1); tt{d}=zeros(sz(d),1);
   return;
end
eps1=eps*nrm/sqrt(d-1);
b=reshape(b,[nl,nr]);
[u0,s0,v0]=svd(b,'econ');
%r=my_chop(diag(s0),eps);
% r=my_chop2(diag(s0),eps1);
r=numel(find(diag(s0)>eps1));
if (exists_max_r) r = min(r, max_r); end;
%r=rank(b,eps1);
u0=u0(:,1:r);
s0=s0(1:r,1:r);
s0=s0./nrm1;
%norm(b-u0*s0*v0(:,1:r)','fro')/norm(b,'fro')
eps1=eps1/nrm1;
b=conj(v0(:,1:r)*s0);
b=conj(b');
tt{1}=u0*nrm1;
for i=2:d-1
    nl=sz(i);
    nr=nr/sz(i);
    b=reshape(b,[r*nl,nr]);
    [u0,s0,v0]=svd(b,'econ');
    %fprintf('%f \n',norm(diag(s0)));
    %r1=my_chop(diag(s0),eps);
    %r1=rank(b,eps1);
%     r1=my_chop2(diag(s0),eps1);
    r1=numel(find(diag(s0)>eps1));
    if (exists_max_r) r1 = min(r1, max_r); end;
    u0=u0(:,1:r1).*nrm1; %u0 is rxnlxr1
    s0=s0(1:r1,1:r1)./nrm1;
    b=conj(v0(:,1:r1)*s0);
    b=conj(b'); %B is r1xn2xn3
    tt{i}=reshape(u0,[r,nl,r1]);
    tt{i}=permute(tt{i},[2,1,3]);   
    r=r1;
    eps1=eps1./nrm1;
end
%In the final, b is just r1xn_d
%fprintf('%f \n',nrm);
tt{d}=conj(b');
return
end
