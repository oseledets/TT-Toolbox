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
    sz=n(1)*rk(2);
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
    sz=n(1)*rk(2);
    cr(1:sz)=cr(1:sz)*a;
    c.core=cr;
    c.r=rk;
    c.n=n;
    c.ps=b.ps;
else
    error('Use mtimes(full(A),full(B)).');
end
