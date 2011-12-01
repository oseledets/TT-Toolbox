function [a]=times(b,c,varargin)
%[A]=TIMES(B,C)
%[A]=TIMES(B,C,EPS)
%Hadamard product of two TT-matrices
%If optional parameter EPS is specified, use Krylov method
if (nargin == 2 )
n = b.n;
m = b.m;
b = b.tt;
c = c.tt;
a = b.*c;
a = tt_matrix(a, n, m);
end
return
end