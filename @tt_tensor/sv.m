function [sval]=sv(tt)
%Computes singular values of all unfoldings of a TT-tensor
%   [SVAL]=SV(TT) Computes all nonzero singular values of all unfoldings 
%    and returns result as a cell array
%    
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

d=tt.d;
n=tt.n;
r=tt.r;
pos=tt.ps;
cr=tt.core;
pos1=1;
sval=cell(d-1,1);
nrm=zeros(d,1);
core0=cr(1:r(1)*n(1)*r(2));
%Orthogonalization from left-to-tight
for i=1:d-1
   core0=reshape(core0,[r(i)*n(i),r(i+1)]);
   [core0,ru]=qr(core0,0); nrm(i+1)=norm(ru,'fro');
   ru=ru./nrm(i+1);
   core1=cr(pos(i+1):pos(i+2)-1);
   core1=reshape(core1,[r(i+1),n(i+1)*r(i+2)]);
   core1=ru*core1;
   r(i+1)=size(core0,2);
   cr(pos1:pos1-1+r(i)*n(i)*r(i+1))=core0(:);
   cr(pos1+r(i)*n(i)*r(i+1):pos1+r(i)*n(i)*r(i+1)+r(i+1)*n(i+1)*r(i+2)-1)=core1(:);
   core0=core1;
   pos1=pos1+r(i)*n(i)*r(i+1);
end
pos1=pos1+r(d)*n(d)*r(d+1)-1;
cr=cr(1:pos1); %Truncate storage if required
pos=cumsum([1;n.*r(1:d).*r(2:d+1)]); 
core0=cr(pos1-r(d)*n(d)*r(d+1)+1:pos1);
 for i=d:-1:2
     %core0=core(pos(i):pos(i+1)-1);
     core1=cr(pos(i-1):pos(i)-1); 
     core0=reshape(core0,[r(i),n(i)*r(i+1)]);
     core1=reshape(core1,[r(i-1)*n(i-1),r(i)]);
     [u,s,v]=svd(core0,'econ');
     s=diag(s); sval{i-1}=s; 
     r(i)=numel(s);
     u=u*diag(s);
     core1=core1*u;
     core0=v';
     cr(pos1-r(i)*n(i)*r(i+1)+1:pos1)=core0(:);
     cr(pos1-r(i)*n(i)*r(i+1)-r(i-1)*n(i-1)*r(i)+1:pos1-r(i)*n(i)*r(i+1))=core1(:);
     %cr=core(pos(i):pos(i+1)-1); 
     pos1=pos1-r(i)*n(i)*r(i+1);
     core0=core1;
 end
 pos1=pos1-r(1)*n(1)*r(2);
 cr=cr(pos1+1:numel(cr)); %Truncate unwanted elements;
 tt.r=r;
 tt.ps=cumsum([1;tt.n.*tt.r(1:d).*tt.r(2:d+1)]);
 pp=cr(1:r(1)*n(1)*r(2));
 nrm(1)=norm(pp,'fro');
 cr(1:r(1)*n(1)*r(2))=pp./nrm(1);
 %Now a simple trick: balance the product of numbers;
 %All cores are orthogonal except the first one. Thus, we know the norm
 nrm0=sum(log(abs(nrm))); 
 nrm0=nrm0/d; nrm0=exp(nrm0);
 %Construct normalization of norm
 for i=1:d-1
   nrm(i+1)=nrm(i+1)*nrm(i)/nrm0;
   nrm(i)=nrm0;
 end
 %Finally redistribute the norm
 ps=tt.ps;
 for i=1:d
    core1=cr(ps(i):ps(i+1)-1);
    core1=core1*nrm(i);
    cr(ps(i):ps(i+1)-1)=core1;
 end
 tt.core=cr;
return
end
