function [y]=mvk(a,x,eps,nswp,z,rmax)
%DMRG fast matrix-by-vector product, Method 1 (less accurate)
%   [Y]=MVK(A,X,EPS,NSWP,Y,RMAX) Matrix-by-vector product of a TT-matrix A
%   by a TT-tensor X with accuracy EPS. Also, one can specify the number of
%   sweeps NSWP, initial approximation Z and the maximal TT-rank RMAX (if
%   they become too large)
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
n=a.n;
m=a.m;
att=a.tt;
corea=att.core;
psa=att.ps;
ra=att.r;
d=att.d;
if ( nargin <= 5 || isempty(rmax) )
   rmax=1000;
end
if ( nargin <= 4 || isempty(z) )
  %z=tt_random(n,ndims(x),5);
  %z=tt_tensor(z);
  z=x;
end
if ( nargin <= 3 || isempty(nswp) )
  nswp = 40;
end
y=z;

%Warmup is to orthogonalize Y from right-to-left and compute psi-matrices
%for Ax
psi=cell(d+1,1); %Psi-matrices 
psi{d+1}=1; psi{1}=1;


%Parameters section

%Warmup: right-to-left QR + computation of psi matrices 

  corex=x.core;
  psx=x.ps;
  rx=x.r;

    psy=y.ps;
  ry=y.r;
  corey=y.core;
   pos1=psy(d+1);
cr1=corey(psy(d):psy(d+1)-1);
for i=d:-1:2  
   cr2=corey(psy(i-1):psy(i)-1);
   cr1=reshape(cr1,[ry(i),n(i)*ry(i+1)]);
   cr2=reshape(cr2,[ry(i-1)*n(i-1),ry(i)]);
   cr1=cr1.';
   [q,rm]=qr(cr1,0); rn=size(q,2); rm=rm.';
   q=q.'; 
   ry(i)=rn;
   corey(pos1-ry(i+1)*n(i)*ry(i):pos1-1)=q(:);
   %Convolution is now performed for psi(i) using psi(i+1) and corea, and
   %(new) core q
   cra=corea(psa(i):psa(i+1)-1); cra=reshape(cra,[ra(i),n(i),m(i),ra(i+1)]);
   cry=reshape(conj(q),[ry(i),n(i),ry(i+1)]);
   crx=corex(psx(i):psx(i+1)-1); crx=reshape(crx,[rx(i),m(i),rx(i+1)]);
   pscur=psi{i+1}; pscur=reshape(pscur,[ra(i+1),rx(i+1),ry(i+1)]); %ra,rx,ry
   %First, convolve over rx(i+1) 
   crx=reshape(crx,[rx(i)*m(i),rx(i+1)]);
   pscur=permute(pscur,[2,1,3]); pscur=reshape(pscur,[rx(i+1),ra(i+1)*ry(i+1)]);
   pscur=crx*pscur; %pscur is now rx(i)*m(i)*ra(i+1)*ry(i+1)
   %Convolve over m(i),ra(i+1),n(i),ry(i+1)
    pscur=reshape(pscur,[rx(i),m(i)*ra(i+1),ry(i+1)]);
    pscur=permute(pscur,[1,3,2]); 
    pscur=reshape(pscur,[rx(i)*ry(i+1),m(i)*ra(i+1)]);
    cra=reshape(cra,[ra(i)*n(i),m(i)*ra(i+1)]); cra=cra.';  
    pscur=pscur*cra; 
    %pscur is now rx(i)*ry(i+1)*ra(i)*n(i), it is left to convolve over 
    %n(i)*ry(i+1)
    pscur=reshape(pscur,[rx(i),ry(i+1),ra(i),n(i)]);
    pscur=permute(pscur,[3,1,4,2]);
    pscur=reshape(pscur,[rx(i)*ra(i),n(i)*ry(i+1)]);
    cry=reshape(cry,[ry(i),n(i)*ry(i+1)]); cry=cry.';
    pscur=pscur*cry;
    psi{i}=pscur;
   %End of psi-block
   pos1=pos1-ry(i+1)*n(i)*ry(i);
   cr1=cr2*rm;
   
end
corey(pos1-ry(2)*n(1)*ry(1):pos1-1)=cr1(:);
pos1=pos1-ry(2)*n(1)*ry(1);
corey=corey(pos1:numel(corey)); %Truncate unused elements

  
psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
y.core=corey;
y.r=ry;
y.ps=psy;
swp=1;
   converged=false;
while (swp <= nswp && ~converged) 
  %left-to-right dmrg sweep 
psy=y.ps;
  ry=y.r;
  corey=y.core;
  
  pos1=1;
  cry_old=corey;
  converged=true;
  ermax=0;
  
  for i=1:d-1
     %We care for two cores, with number i & number i+1, and use
     %psi(i) and psi(i+2) as a basis; also we will need to recompute
     %psi(i+1)
     ps1=psi{i}; ps2=psi{i+2};
     cra1=corea(psa(i):psa(i+1)-1); cra2=corea(psa(i+1):psa(i+2)-1);
     crx1=corex(psx(i):psx(i+1)-1); crx2=corex(psx(i+1):psx(i+2)-1);
     %our convolution is
     %ps1(ra(i),rx(i),ry(i))*cra1(ra(i),n(i),m(i),ra(i+1))*
     %cra2(ra(i+1),n(i+1),m(i+1),ra(i+2))*
     %*ps2(ra(i+2),rx(i+2),ry(i+2))
     %*crx1(rx(i),m(i),rx(i+1))*cr2x(rx(i+1)*m(i+1)*rx(i+2))
     %Scheme ps1*crx1 over rx(i)
     
     ps1=reshape(ps1,[ra(i),rx(i),ry(i)]); 
     ps1=permute(ps1,[1,3,2]); ps1=reshape(ps1,[ra(i)*ry(i),rx(i)]);
     crx1=reshape(crx1,[rx(i),m(i)*rx(i+1)]);
     ps1=ps1*crx1; %ps1 is now ra(i)*ry(i)*m(i)*rx(i+1)
     %Now convolve with matrix A over ra(i)*m(i)
     ps1=reshape(ps1,[ra(i),ry(i),m(i),rx(i+1)]);
     ps1=permute(ps1,[2,4,1,3]);
     ps1=reshape(ps1,[ry(i)*rx(i+1),ra(i)*m(i)]);
     cra1=reshape(cra1,[ra(i),n(i),m(i),ra(i+1)]);
     cra1=permute(cra1,[1,3,2,4]); 
     cra1=reshape(cra1,[ra(i)*m(i),n(i)*ra(i+1)]);
     ps1=ps1*cra1; %ps1 is now ry(i)*rx(i+1)*n(i)*ra(i+1)
     %Then the ``same'' convolution is carried over for second pair of
     %cores
     ps2=reshape(ps2,[ra(i+2),rx(i+2),ry(i+2)]);
     ps2=permute(ps2,[2,1,3]);
     ps2=reshape(ps2,[rx(i+2),ra(i+2)*ry(i+2)]);
     crx2=reshape(crx2,[rx(i+1)*m(i+1),rx(i+2)]);
     ps2=crx2*ps2; %ps2 is now rx*(i+1)*m(i+1)*ra(i+2)*ry(i+2)
     %Convolve over m(i+1)*ra(i+2)
     ps2=reshape(ps2,[rx(i+1),m(i+1)*ra(i+2),ry(i+2)]);
     ps2=permute(ps2,[2,1,3]);
     ps2=reshape(ps2,[m(i+1)*ra(i+2),rx(i+1)*ry(i+2)]);
     cra2=reshape(cra2,[ra(i+1)*n(i+1),m(i+1)*ra(i+2)]);
     ps2=cra2*ps2; 
     %ps2 is now ra(i+1)*n(i+1)*rx(i+1)*ry(i+2)
     %Now form superblock by contraction ps1 & ps2
     %over ra(i+1)*rx(i+1)
     ps0=reshape(ps1,[ry(i),rx(i+1),n(i),ra(i+1)]);
     ps0=permute(ps0,[1,3,2,4]);
     ps0=reshape(ps0,[ry(i)*n(i),rx(i+1)*ra(i+1)]);
     ps2=reshape(ps2,[ra(i+1),n(i+1),rx(i+1),ry(i+2)]);
     ps2=permute(ps2,[3,1,2,4]); 
     ps2=reshape(ps2,[rx(i+1)*ra(i+1),n(i+1)*ry(i+2)]);
     super_core=ps0*ps2; %super_core is ry(i)*n(i)*n(i+1)*ry(i+2)
     %Compute previous supercore
     cr1=corey(pos1:pos1+ry(i)*n(i)*ry(i+1)-1);
     cr2=cry_old(psy(i+1):psy(i+2)-1); 
     cr1=reshape(cr1,[ry(i)*n(i),ry(i+1)]);
     cr2=reshape(cr2,[ry(i+1),n(i+1)*ry(i+2)]);
     super_core_old=cr1*cr2; 
     er=norm(super_core_old(:)-super_core(:))/norm(super_core(:));
     if ( er > eps) 
        converged=false;
     end
     ermax=max(er,ermax);
     %fprintf('i=%d er=%3.2e \n',i,er);
     [u,s,v]=svd(super_core,'econ');
     s=diag(s); 
     r=my_chop2(s,eps/sqrt(d-1)*norm(s)); r=min(r,rmax);
     u=u(:,1:r); s=s(1:r); v=v(:,1:r); v=v*diag(s);
     ry(i+1)=r;
     %u is ry(i)*n(i)*ry(i+1)
     %core_new(pos1:pos1+ry(i)*n(i)*ry(i+1)-1)=u(:); 
     corey(pos1:pos1+ry(i)*n(i)*ry(i+1)-1)=u(:); 
     
     u=reshape(u,[ry(i),n(i),ry(i+1)]);

     %Compute new psi
     %ps1 is ry(i)*rx(i+1)*n(i)*ra(i+1) with u over ry(i)*n(i)
     u=conj(u);
     ps1=reshape(ps1,[ry(i),rx(i+1),n(i),ra(i+1)]);
     ps1=permute(ps1,[4,2,1,3]);
     ps1=reshape(ps1,[ra(i+1)*rx(i+1),ry(i)*n(i)]);
     u=reshape(u,[ry(i)*n(i),ry(i+1)]);
     ps1=ps1*u;
     psi{i+1}=ps1;
     %Compute (?) new v
     pos1=pos1+ry(i)*n(i)*ry(i+1);
     v=v';
     %v=reshape(v,[ry(i+1),n(i+1),ry(i+2)]);
     corey(pos1:pos1+ry(i+1)*n(i+1)*ry(i+2)-1)=v(:);
  end
  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
y.core=corey;
y.r=ry;
y.ps=psy;
cry_old=corey;
%m=numel(cry);
%cry=[zeros(size(cry)),cry];
%pos1=pos1+m;
%cry2=cry_old(psy(i+1):psy(i+2)-1);
corey=corey(psy(d):psy(d+1)-1); %Start--only two cores are left!

 
  for i=d-1:-1:1
     %We care for two cores, with number i & number i+1, and use
     %psi(i) and psi(i+2) as a basis; also we will need to recompute
     %psi(i)
     ps1=psi{i}; ps2=psi{i+2};
     cra1=corea(psa(i):psa(i+1)-1); cra2=corea(psa(i+1):psa(i+2)-1);
     crx1=corex(psx(i):psx(i+1)-1); crx2=corex(psx(i+1):psx(i+2)-1);
     %our convolution is
     %ps1(ra(i),rx(i),ry(i))*cra1(ra(i),n(i),m(i),ra(i+1))*
     %cra2(ra(i+1),n(i+1),m(i+1),ra(i+2))*
     %*ps2(ra(i+2),rx(i+2),ry(i+2))
     %*crx1(rx(i),m(i),rx(i+1))*cr2x(rx(i+1)*m(i+1)*rx(i+2))
     %Scheme ps1*crx1 over rx(i)
     
     ps1=reshape(ps1,[ra(i),rx(i),ry(i)]); 
     ps1=permute(ps1,[1,3,2]); ps1=reshape(ps1,[ra(i)*ry(i),rx(i)]);
     crx1=reshape(crx1,[rx(i),m(i)*rx(i+1)]);
     ps1=ps1*crx1; %ps1 is now ra(i)*ry(i)*m(i)*rx(i+1)
     %Now convolve with matrix A over ra(i)*m(i)
     ps1=reshape(ps1,[ra(i),ry(i),m(i),rx(i+1)]);
     ps1=permute(ps1,[2,4,1,3]);
     ps1=reshape(ps1,[ry(i)*rx(i+1),ra(i)*m(i)]);
     cra1=reshape(cra1,[ra(i),n(i),m(i),ra(i+1)]);
     cra1=permute(cra1,[1,3,2,4]); 
     cra1=reshape(cra1,[ra(i)*m(i),n(i)*ra(i+1)]);
     ps1=ps1*cra1; %ps1 is now ry(i)*rx(i+1)*n(i)*ra(i+1)
     %Then the ``same'' convolution is carried over for second pair of
     %cores
     ps2=reshape(ps2,[ra(i+2),rx(i+2),ry(i+2)]);
     ps2=permute(ps2,[2,1,3]);
     ps2=reshape(ps2,[rx(i+2),ra(i+2)*ry(i+2)]);
     crx2=reshape(crx2,[rx(i+1)*m(i+1),rx(i+2)]);
     ps2=crx2*ps2; %ps2 is now rx*(i+1)*m(i+1)*ra(i+2)*ry(i+2)
     %Convolve over m(i+1)*ra(i+2)
     ps2=reshape(ps2,[rx(i+1),m(i+1)*ra(i+2),ry(i+2)]);
     ps2=permute(ps2,[2,1,3]);
     ps2=reshape(ps2,[m(i+1)*ra(i+2),rx(i+1)*ry(i+2)]);
     cra2=reshape(cra2,[ra(i+1)*n(i+1),m(i+1)*ra(i+2)]);
     ps2=cra2*ps2; 
     %ps2 is now ra(i+1)*n(i+1)*rx(i+1)*ry(i+2)
     %Now form superblock by contraction ps1 & ps2
     %over ra(i+1)*rx(i+1)
     ps0=reshape(ps1,[ry(i),rx(i+1),n(i),ra(i+1)]);
     ps0=permute(ps0,[1,3,2,4]);
     ps0=reshape(ps0,[ry(i)*n(i),rx(i+1)*ra(i+1)]);
     ps2=reshape(ps2,[ra(i+1),n(i+1),rx(i+1),ry(i+2)]);
     ps2=permute(ps2,[3,1,2,4]); 
     ps2=reshape(ps2,[rx(i+1)*ra(i+1),n(i+1)*ry(i+2)]);
     super_core=ps0*ps2; %super_core is ry(i)*n(i)*n(i+1)*ry(i+2)
     %Compute previous supercore
     cry1=cry_old(psy(i):psy(i+1)-1);
     %cry2=cry(ry(:ry(i+1)*n(i+1)*ry(i+2));
     %cry1=cry(1:ry(i)*n(i)*ry(i+1));
     %cry2=cry(ry(i)*n(i)*ry(i+1)+1:ry(i)*n(i)*ry(i+1)+ry(i+1)*n(i+1)*ry(i+2));
     cry2=corey(1:ry(i+1)*n(i+1)*ry(i+2));
     corey(1:ry(i+1)*n(i+1)*ry(i+2))=[];
     cry1=reshape(cry1,[ry(i)*n(i),ry(i+1)]);
     cry2=reshape(cry2,[ry(i+1),n(i+1)*ry(i+2)]);
     super_core_old=cry1*cry2;
     er=norm(super_core_old(:)-super_core(:))/norm(super_core(:));
     if ( er > eps) 
        converged=false;
     end
     er=max(er,ermax);
     [u,s,v]=svd(super_core,'econ');
     s=diag(s); 
     r=my_chop2(s,eps/sqrt(d-1)*norm(s)); r=min(r,rmax);
     u=u(:,1:r); s=s(1:r); u=u*diag(s); v=v(:,1:r);
     ry(i+1)=r;
     v=v.';
       u1=u(:); u1=u1.'; v1=v(:); v1=v1.';
     corey=[u1,v1,corey];
      %Compute new psi
      %ps2 is rx(i+1)*ra(i+1)*n(i+1)*ry(i+2)
      %v is ry(i+1)*n(i+1)*ry(i+2)
      %convolve over n(i+1)*ry(i+2);
      v = v.';
      ps2=reshape(ps2,rx(i+1)*ra(i+1),n(i+1)*ry(i+2));
      ps2=ps2*v; %ps2 is now rx(i+1)*ra(i+1)*ry(i+1)
      ps2=reshape(ps2,[rx(i+1),ra(i+1),ry(i+1)]);
      ps2=permute(ps2,[2,1,3]);
      psi{i+1}=ps2;
  end  
  swp=swp+1;
  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
  %fprintf('swp=%d er=%3.2e \n',swp,ermax);
  y.core=corey;
  y.r=ry;
  y.ps=psy;
end
%Compute alpha*y0 = Ax - y
%p=p/tt_dot(y0,y0);
%keyboard;%
%end
%p1=tt_mvdot(core(a),core(x),y0);
%p2=dot(y,tt_tensor(y0));
%abs(p1-p2)/abs(p1)
%keyboard;


  