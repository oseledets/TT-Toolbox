function [tt]=round(tt,varargin)
%Approximate TT-tensor with another one with specified accuracy
%   [TT]=ROUND(TT,EPS) Approximate TT-tensor with relative accuracy EPS
%
%   [TT]=ROUND(TT,EPS,RMAX) Approximate TT-tensor with relative accuracy 
%   EPS and maximal rank RMAX. RMAX can be array of ranks or a number
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

%
if (nargin == 2 )
  eps=varargin{1};
  rmax=prod(size(tt));
elseif ( nargin == 3 )
  eps=varargin{1};
  rmax=varargin{2};
end
d=tt.d;
n=tt.n;
r=tt.r;
if (numel(rmax) == 1 )
  rmax=rmax*ones(d+1,1); 
end
pos=tt.ps;
cr=tt.core;
pos1=1;
nrm=zeros(d,1);
core0=cr(1:r(1)*n(1)*r(2));
%Orthogonalization from left-to-tight
for i=1:d-1
   core0=reshape(core0,[r(i)*n(i),r(i+1)]);
   [core0,ru]=qr(core0,0); nrm(i+1)=norm(ru,'fro');
   if (nrm(i+1)~=0)
    ru=ru./nrm(i+1);
   end;
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
ep=eps/sqrt(d-1);
pos=cumsum([1;n.*r(1:d).*r(2:d+1)]); 
core0=cr(pos1-r(d)*n(d)*r(d+1)+1:pos1);
for i=d:-1:2
     core1=cr(pos(i-1):pos(i)-1); 
     core0=reshape(core0,[r(i),n(i)*r(i+1)]);
     core1=reshape(core1,[r(i-1)*n(i-1),r(i)]);

     try
       [u,s,v]=svd(core0, 'econ');
     catch err
       if (strcmp(err.identifier,'MATLAB:svd:svdNoConvergence'))
         % When svd doesn't converge, svds usually works fine (but much slower).
         warning('SVD did not converge, switched to svds.');
         [u,s,v]=svds(core0, min(size(core0)));
       else
         rethrow(err);
       end
     end

     s=diag(s);
     r1=my_chop2(s,norm(s)*ep);
     r1=min(r1,rmax(i));

     u=u(:,1:r1);s=s(1:r1); v=v(:,1:r1);
     u=u*diag(s);
     r(i)=r1;
     core1=core1*u;
     core0=v';
     cr(pos1-r(i)*n(i)*r(i+1)+1:pos1)=core0(:);
     cr(pos1-r(i)*n(i)*r(i+1)-r(i-1)*n(i-1)*r(i)+1:pos1-r(i)*n(i)*r(i+1))=core1(:);
     pos1=pos1-r(i)*n(i)*r(i+1);
     core0=core1;
end
pos1 = pos1-r(1)*n(1)*r(2);
cr = cr(pos1+1:numel(cr)); %Truncate unwanted elements;
tt.r=r;
tt.ps=cumsum([1;tt.n.*tt.r(1:d).*tt.r(2:d+1)]);
pp=cr(1:r(1)*n(1)*r(2));
nrm(1) = norm(pp,'fro');
if (nrm(1) ~= 0)
    pp = pp./nrm(1);
end;
cr(1:r(1)*n(1)*r(2)) = pp;
 %Now a simple trick: balance the product of numbers;
 %All cores are orthogonal except the first one. Thus, we know the norm
 nrm0=sum(log(abs(nrm))); 
 nrm0=nrm0/d; nrm0=exp(nrm0);
 ps=tt.ps;
 for i=1:d
    cr(ps(i):ps(i+1)-1)=nrm0*cr(ps(i):ps(i+1)-1);
 end
 tt.core=cr;
 tt.over=0;
return
end
