function [c] = mtimes(a,b)
%[C]=MTIMES(A,B)
%Implement multiplication for a tt_tensor and (in future) 
%matrix-by-vector product

if isa(a,'tt_tensor') && numel(b) == 1
    c = tt_tensor; 
    c.d = a.d;
    n = a.n;
    rk=a.r;
    cr=a.core;
    sz=rk(1)*n(1)*rk(2);
    cr(1:sz)=cr(1:sz)*b;
    c.core=cr;
    c.r=rk;
    c.n=n;
    c.ps=a.ps;
elseif isa(b,'tt_tensor') && numel(a) == 1
    c = tt_tensor; 
    c.d = b.d;
    n = b.n;
    rk=b.r;
    cr=b.core;
    sz=rk(1)*n(1)*rk(2);
    cr(1:sz)=cr(1:sz)*a;
    c.core=cr;
    c.r=rk;
    c.n=n;
    c.ps=b.ps;
elseif isa(a,'double') && isa(b,'tt_tensor')
    c=b;
    crc=c.core; 
    if (size(crc,2) == 1 ) 
      crc=crc.';
    end
    ps=c.ps;
    r=c.r;
    n=c.n;
    mm=crc(ps(1):ps(2)-1);
    mm=reshape(mm,[r(1),n(1)*r(2)]);
    mm=a*mm;
    crc(ps(1):ps(2)-1)=[];
    crc=[mm(:)',crc];
    r(1)=size(a,1);
    ps=ps-ps(2)+1+r(1)*n(2)*r(2);
    ps(1)=1;
    c.ps=ps;
    c.core=crc;
    c.r=r;
else
    error('Use mtimes(full(A),full(B)).');
end
