function [c] = kron(a,b)
%Kronecker product of two TT-tensors, canonical type
%   [C]=KRON(A,B) Kronecker product of two TT-tensors. One of the
%   arguments can be empty. Result has the same dimension as the arguments.
%   Mode sizes are multiplied.
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

if isempty(a)
    c=b;
    return;
end
if isempty(b)
    c=a;
    return;
end
if ndims(a) == ndims(b)
    c = tt_tensor;
    na = a.n; nb = b.n;
    c.n=na.*nb;
    d = ndims(a);
    ra = a.r; rb = b.r;
    c.d = d;
    c.r = ra.*rb;
    c.ps = cumsum([1; c.r(1:d).*c.n.*c.r(2:d+1)]);
    c.core = zeros(c.ps(d+1)-1, 1);
    for i=1:d
        cra = reshape(a.core(a.ps(i):a.ps(i+1)-1), 1, ra(i)*na(i)*ra(i+1));
        crb = reshape(b.core(b.ps(i):b.ps(i+1)-1), rb(i)*nb(i)*rb(i+1), 1);
        crc = crb*cra; % Canonical big-endian
        crc = reshape(crc, rb(i), nb(i), rb(i+1), ra(i), na(i), ra(i+1));
        crc = permute(crc, [1, 4, 2, 5, 3, 6]);
        c.core(c.ps(i):c.ps(i+1)-1) = crc(:);
    end
else
    error('Tensor dimensions must agree');
end

end
