function [c] = kron(a,b)
%Kronecker product of two TT-matrices, canonical type
%   [C]=KRON(A,B) Kronecker product of two TT-matrices. One of the
%   arguments can be empty. Result has the same dimension as the arguments.
%   Mode sizes are multiplied.
%
% TT-Toolbox 2.2, 2009-2013
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

if isempty(a)
    c=b;
    return;
end
if isempty(b)
    c=a;
    return;
end
if ndims(a.tt) == ndims(b.tt)
    c=tt_matrix;
    na = a.n; nb = b.n;
    ma = a.m; mb = b.m;
    c.n=na.*nb;
    c.m=ma.*mb; 
    d = ndims(a.tt);
    ra = a.tt.r; rb = b.tt.r;
    c.tt.d = d;
    c.tt.r = ra.*rb;
    c.tt.n = c.m.*c.n;
    c.tt.ps = cumsum([1; c.tt.r(1:d).*c.tt.n.*c.tt.r(2:d+1)]);
    c.tt.core = zeros(c.tt.ps(d+1)-1, 1);
    for i=1:d
        cra = reshape(a.tt{i}, [ra(i)*na(i),ma(i)*ra(i+1)]);
        crb = reshape(b.tt{i}, [rb(i)*nb(i),mb(i)*rb(i+1)]);
        crc = kron(crb,cra);
        crc = reshape(crc, [ra(i),na(i),rb(i),nb(i),ma(i),ra(i+1),mb(i),rb(i+1)]);
        crc = permute(crc, [1, 3, 4, 2, 7, 5, 6, 8]);
        c.tt.core(c.tt.ps(i):c.tt.ps(i+1)-1) = crc(:);
    end    
else
    error('Tensor dimensions must agree');
end
return
end
