function [y]=axpx(a,x)
%[y]=AXPX(A,X)
%This is a technical function that multiplies
%Matrix A by a vector and projects on X itself; 
%This is for time-stepping 

d=x.d;
n=a.n;
m=a.m;
rx=x.r;
psx=x.ps;
att=a.tt;
corea=att.core;
psa=att.ps;
ra=att.r;
%Do right-to-left orthogonalization and computation of psi
%matrices
corex=x.core;
pos1=psx(d+1);
psi=cell(d+1,1); %Psi-matrices 
psi{d+1}=1; psi{1}=1;
cr1=corex(psx(d):psx(d+1)-1);
y=x;
corey=y.core;
psy=y.ps;
ry=y.r;
for i=d:-1:2 
   %  fprintf('i=%d \n');
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
   crx=corex(psx(i):psx(i+1)-1);
   crx=reshape(crx,[rx(i),m(i),rx(i+1)]);
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
    pscur=reshape(pscur,[ra(i)*rx(i),n(i)*ry(i+1)]);
    cry=reshape(cry,[ry(i),n(i)*ry(i+1)]); cry=cry.';
    pscur=pscur*cry; %pscur is ra*rx*ry
    pscur=reshape(pscur,[rx(i),ra(i),ry(i)]);
    psi{i}=pscur;
   %End of psi-block
   pos1=pos1-ry(i+1)*n(i)*ry(i);
   cr1=cr2*rm;
end
corey(pos1-ry(2)*n(1)*ry(1):pos1-1)=cr1(:);
pos1=pos1-ry(2)*n(1)*ry(1);
corey=corey(pos1:end); %Truncate unused elements
%Now compute the projection itself while recomputing psi-matrices
 pos1=1;
 psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
 rm=1;
 %keyboard;
  for i=1:d
          %Our convolution is
     %ps1(ra(i),rx(i),ry(i))*cra(ra(i),n(i),m(i),ra(i+1))
     %*ps2(ra(i+1),rx(i+1),ry(i+1))*crx(rx(i),m(i)*rx(i+1))->
     %cry(ry(i),n(i),ry(i+1)

     ps1=psi{i}; ps2=psi{i+1};
     cra=corea(psa(i):psa(i+1)-1);
     crx=corex(psx(i):psx(i+1)-1);
%     corey(psy(i):psy(i+1)-1)=crx(:);
%     y.core=corey;
%          psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
%     y.ps=psy;
     ps1=reshape(ps1,[ra(i),rx(i),ry(i)]);
     ps1=permute(ps1,[2,3,1]);
     ps1=reshape(ps1,[rx(i)*ry(i),ra(i)]);
     cra=reshape(cra,[ra(i),n(i)*m(i)*ra(i+1)]);
     cry=ps1*cra;      %cry is rx(i)*ry(i)*n(i)*m(i)*ra(i+1)
     %convolve over m(i),rx(i) with crx
     cry=reshape(cry,[rx(i),ry(i),n(i),m(i),ra(i+1)]);
     cry=permute(cry,[2,3,5,1,4]);
     cry=reshape(cry,[ry(i)*n(i)*ra(i+1),rx(i)*m(i)]);
     crx=reshape(crx,[rx(i)*m(i),rx(i+1)]);
     cry=cry*crx; cry0=cry;
     %cry is ry(i)*n(i)*ra(i+1)*rx(i+1)
     cry=reshape(cry,[ry(i)*n(i),ra(i+1)*rx(i+1)]);
     ps2=reshape(ps2,[ra(i+1)*rx(i+1),ry(i+1)]);
     %cry0=cry; %cry0 is ry(i)*n(i)*ra(i+1)*rx(i+1), convolution over ry(i)*n(i) will bring new psi --- not true; new psi will be brought with and
     %orthogonalization
     cry=cry*ps2;
      corey(pos1:pos1+ry(i)*n(i)*ry(i+1)-1)=cry(:);
      y.core=corey;
      psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
      y.ps=psy;
       %keyboard;
    
     cry=reshape(cry,[ry(i)*n(i),ry(i+1)]);
     [q,rm]=qr(cry,0);
     rn=size(q,2);
     ry(i+1)=rn;
     %q=q*rm;
     corey(pos1:pos1+ry(i)*n(i)*ry(i+1)-1)=q(:);
     pos1=pos1+ry(i)*n(i)*ry(i+1);
     %And now --- orthogonalization! (of cry) (he-he)
     %Now compute "new" psi
     %We need to compute 
          %cry is rx(i)*ry(i)*n(i)*m(i)*ra(i+1)
     %convolve over m(i),rx(i) with crx
     %cry0 is ry(i)*n(i)*ra(i+1)*rx(i+1); 
     cry0=reshape(cry0,[ry(i),n(i),ra(i+1),rx(i+1)]);
     cry0=permute(cry0,[3,4,1,2]);
     cry0=reshape(cry0,[ra(i+1)*rx(i+1),ry(i)*n(i)]);
     q=reshape(conj(q),[ry(i)*n(i),ry(i+1)]);
     psi{i+1}=cry0*q;
     
  end
  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
  corey(psy(i):psy(i+1)-1)=corey(psy(i):psy(i+1)-1)*rm;
  y.core=corey;
  y.r=ry;
  y.ps=psy;
  y.n=n;
  y.d=d;
  return
end