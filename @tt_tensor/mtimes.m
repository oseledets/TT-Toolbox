function [c] = mtimes(a,b)
%   C = A*B
%   [C]=MTIMES(A,B)
%   TT-Tensor-by-Matrix multiplication
%
%   Either A or B may be a TT-tensor. The rest argument must be a matrix
%   with sizes consistent to the respective "tail" rank of a tensor (or a
%   scalar).

if isa(a,'double') && isa(b,'tt_tensor')
    d = b.d;
    cr1 = b.core(b.ps(1):b.ps(2)-1);
    cr1 = reshape(cr1, b.r(1), b.n(1)*b.r(2));
    if (numel(a)==1) 
        rnew = b.r(1);
    else
        rnew = size(a,1);
    end;
    dps = rnew*b.n(1)*b.r(2) - b.r(1)*b.n(1)*b.r(2);
    cr1 = a*cr1;
    
    c = b;
    c.core=c.core(b.ps(2):b.ps(d+1)-1);
    c.core = [cr1(:); c.core(:)];
    c.ps(2:d+1)=b.ps(2:d+1)+dps;
    c.r(1) = rnew;
    
%     c=b;
%     crc=c.core; 
%     if (size(crc,2) == 1 ) 
%       crc=crc.';
%     end
%     ps=c.ps;
%     r=c.r;
%     n=c.n;
%     mm=crc(ps(1):ps(2)-1);
%     mm=reshape(mm,[r(1),n(1)*r(2)]);
%     mm=a*mm;
%     crc(ps(1):ps(2)-1)=[];
%     crc=[mm(:)',crc];
%     r(1)=size(a,1);
%     ps=ps-ps(2)+1+r(1)*n(2)*r(2);
%     ps(1)=1;
%     c.ps=ps;
%     c.core=crc;
%     c.r=r;
elseif isa(a,'tt_tensor') && isa(b,'double') 
    d = a.d;
    crd = a.core(a.ps(d):a.ps(d+1)-1);
    crd = reshape(crd, a.r(d)*a.n(d), a.r(d+1));
    if (numel(b)==1) 
        rnew = a.r(d+1);
    else
        rnew = size(b,2);
    end;    
    ps_new = a.ps(d) + a.r(d)*a.n(d)*rnew;
    crd = crd*b;
    
    c = a;
    c.core(a.ps(d):ps_new-1)=crd(:);
    c.core = c.core(1:ps_new-1);
    c.ps(d+1)=ps_new;
    c.r(d+1) = rnew;
    
else
    error('Use mtimes(full(A),full(B)).');
end
