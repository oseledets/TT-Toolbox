function [c] = kron(a,b)
%[C]=KRON(A,B)
%Kronecker product of two TT-matrices
c=tt_matrix;
c.n=[(a.n);(b.n)];
c.m=[(a.m);(b.m)];
%c.d=a.d+b.d;
c.tt=kron(a.tt,b.tt);
return
end
