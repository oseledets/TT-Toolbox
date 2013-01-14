function [trA] = trace(a)
%Kronecker product of two TT-matrices
%   [trA]=TRACE(A) Trace of TT-matrix, which should have square dimensions.
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

trA=0.e0;
if isempty(a)
    return;
end
d=ndims(a);
if a.n == a.m
    b=tt_tensor;
    b.d=d;
    b.r=a.tt.r;
    b.n=ones(d,1);
    b.ps = cumsum([1; b.r(1:d).*b.n.*b.r(2:d+1)]);
    b.core = zeros(b.ps(d+1)-1, 1);
    for i=1:d
        cra = reshape(a.tt{i}, [a.tt.r(i),a.n(i)*a.m(i),a.tt.r(i+1)]);
        cra = permute(cra, [2,1,3]);
        cra = reshape(cra, [a.n(i),a.m(i),a.tt.r(i)*a.tt.r(i+1)]);
        crb = zeros(a.tt.r(i)*a.tt.r(i+1),1);
        parfor j=1:a.tt.r(i)*a.tt.r(i+1)
            crb(j)=trace(cra(:,:,j));
        end
        b.core(b.ps(i):b.ps(i+1)-1) = crb(:);
    end
    trA=full(b);
else
    error('Matrix dimensions of A must agree');
end
return
end
