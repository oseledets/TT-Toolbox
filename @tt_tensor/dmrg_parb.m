function [y]=dmrg_parb(mat,a,f,eps,y0,rmax,nswp)
%[y]=dmrg_parb(mat,a,rhs,eps,[y0],[rmax],[nswp])
%This is not an east stuff.
%It solve parameter-dependent problem A(y)u(y) = rhs(y)
%With accuracy eps.
%A(y) is given as a sum mat{s}*a(s,P), and a is stored
%as a row TT-tensor (i.e. r(1) = r
%mat is stored as 1D TT-matrix with r(2) = r

%The storage is a delicate issue. We will stick to the
%separate storage of the "Physical" part and of the diagonal part
%We will focus on the "maxvol" part.


%This is the "physical dimension" size
N=mat.n;
n=a.n;
m=a.m; %It maybe useful
d=a.d; %And this is for the parametric part
%We seek for the solution as JxY1xQ1xQ1xY2xQ2
k=mat.r(2); %This is "the block size"
if ( nargin < 4 || isempty(y0) )
   y=tt_tensor;
   y.d=d;
   ry=k*ones(1,d); ry(1)=k; ry(d+1)=1; ry=ry';
   y.r=ry;
   y.n=n;
   y.ps=cumsum([1;y.n.*ry(1:d).*ry(2:d+1)]);
   psy=y.ps;
   sz=psy(d+1)-1;
   cr=randn(1,sz);
   %cr=ones(sz,1);
   y.core=cr;
else 
   y=y0;
end
if ( nargin < 5 || isempty(rmax) )
  rmax=1000;
end
if ( nargin < 6 || isempty(nswp) )
  nswp=10;
end
%Matrix details
cra=a.core;
ra=a.r;
psa=a.ps;
ry=y.r;
psy=y.ps;
crf=f.core;
psf=f.ps;
rf=f.r;

%Right hand side details --- we would like to split 
%the physical part and the parametric part here also,
%but we do not want to bother the reader with such amount
%of useless parameters, so we will do it ourselves(?)
%NO! We can just compute phi up to i=2 and that is all




%The situation is "kinda" similar; however, we start from right-to-left 
%qr & maxvol of Y; The "matrix" core is like ra(i)xn(i)xm(i)*ra(i+1)
%That means, that the "left" matrix is ry(i)xkxk, and the right is just
%ra(i+1)
%Actually, we need for a specific index subset (y1,...,yk,J,yk+1...yP)
%Compute the matrix; it is more (or less) reduces to the computation
%of the coefficient; i.e. B(i,j,Q)*Phi(Q,ALLINDICES) -> B(I,J,ALLINDICES)
%B(i,j,Q1)*C(Q1,y1,Q2)*C(Q2,y2,Q3)*C(Q3,y3,Q4)*C(Q4,y4) 
%B(i,j,Q1)*P(Q1,Q3,RY1)*C(Q3,y3,Q4)*P(Q4,RY2,Q5)
%B(i,j,Q3,RY1)
rm=1;
pos1=psy(d+1)-1;
phi=cell(d+1,1);
phi{d+1}=1;
phif=cell(d+1,1);
phif{d+1}=1;
%Right-to-left qr & maxvol
for i=d:-1:1
    cr=cry(psy(i):psy(i+1)-1);
    cr=reshape(cr,[ry(i)*n(i),ry(i+1)]);
    cr=cr*rm; 
    ry(i+1)=size(cr,2); 
    cr=reshape(cr,[ry(i),n(i)*ry(i+1)]);
    cr=cr.'; %It seems to be needed
    [u,rm]=qr(cr,0); 
    indu=maxvol2(u); 
    r1=u(indu,:);
    u=u/r1;
    rm=r1*rm;
    u=u.'; %n(i)*ry(i+1)xry(i)->ry(i)xn(i)xry(i+1)
    rnew=size(u,2); 
    %ry(i+1)=rnew;
    
    cry(pos1-ry(i)*n(i)*rnew:pos1)=u(:); 
    pos1=pos1-ry(i)*n(i)*rnew; %Mda
    %With this new core compute new phi matrix; we better need raxry
    
    crm=cra(psa(i):psa(i+1)-1);
    phx=phi{i+1};
    phx=reshape(phx,[ra(i+1),ry(i+1)]);
    crm=reshape(crm,[ra(i)*n(i),ra(i+1)]);
    crm=crm*phx; %crm is ra(i)xn(i)xry(i+1); u is rnewxn(i)xry(i+1);
    crm=reshape(crm,[ra(i),n(i)*ry(i+1)]); 
    phi{i}=crm(:,indu);
    crh=crf(psf(i+1):psi(i+2)-1); 
    crh=reshape(crh,[rf(i+1)*n(i),rf(i+2)]);
    phf=reshape(phif{i+1},[rf(i+2),ry(i+1)]);
    crh=crh*phf;
    crh=reshape(crh,[rf(i+1),n(i)*ry(i+1)]);
    phif{i}=crh(:,indu);
end

%And now truncate pos
cry=cry(pos1:numel(cry));
y.core=cry;
y.r=ry;
y.ps=psy;
%keyboard
swp=1;



return
end