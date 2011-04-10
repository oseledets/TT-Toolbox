function [tt]=round(tt,varargin)
%[TT]=ROUND(TT,EPS)
%Approximate TT-tensor with relative accuracy EPS
%[TT]=ROUND(TT,EPS,RMAX)
%Approximate TT-tensor with relative accuracy EPS and maximal rank rmax
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
pos=tt.ps;
core=tt.core;
pos1=1;
core0=core(1:r(1)*n(1)*r(2));
%Orthogonalization from left-to-tight
for i=1:d-1
   core0=reshape(core0,[r(i)*n(i),r(i+1)]);
   [core0,ru]=qr(core0,0);
   core1=core(pos(i+1):pos(i+2)-1);
   core1=reshape(core1,[r(i+1),n(i+1)*r(i+2)]);
   core1=ru*core1;
   r(i+1)=size(core0,2);
   core(pos1:pos1-1+r(i)*n(i)*r(i+1))=core0(:);
   core(pos1+r(i)*n(i)*r(i+1):pos1+r(i)*n(i)*r(i+1)+r(i+1)*n(i+1)*r(i+2)-1)=core1(:);
   core0=core1;
   pos1=pos1+r(i)*n(i)*r(i+1);
end

pos1=pos1+r(d)*n(d)*r(d+1)-1;
core=core(1:pos1); %Truncate storage if required
 ep=eps/sqrt(d-1);
pos=cumsum([1;n.*r(1:d).*r(2:d+1)]); 
core0=core(pos1-r(d)*n(d)*r(d+1)+1:pos1);
 for i=d:-1:2
     %core0=core(pos(i):pos(i+1)-1);
     core1=core(pos(i-1):pos(i)-1); 
     core0=reshape(core0,[r(i),n(i)*r(i+1)]);
     core1=reshape(core1,[r(i-1)*n(i-1),r(i)]);
     [u,s,v]=svd(core0,'econ');
     s=diag(s); r1=my_chop2(s,norm(s)*ep);
     r1=min(r1,rmax);
     u=u(:,1:r1);s=s(1:r1); v=v(:,1:r1);
     u=u*diag(s);
     r(i)=r1;
     core1=core1*u;
     core0=v';
     core(pos1-r(i)*n(i)*r(i+1)+1:pos1)=core0(:);
     core(pos1-r(i)*n(i)*r(i+1)-r(i-1)*n(i-1)*r(i)+1:pos1-r(i)*n(i)*r(i+1))=core1(:);
     %cr=core(pos(i):pos(i+1)-1); 
     pos1=pos1-r(i)*n(i)*r(i+1);
     core0=core1;
 end
 pos1=pos1-r(1)*n(1)*r(2);
 core=core(pos1+1:numel(core)); %Truncate unwanted elements
 tt.core=core;
 tt.r=r;
 tt.ps=cumsum([1;tt.n.*tt.r(1:d).*tt.r(2:d+1)]);
return
end