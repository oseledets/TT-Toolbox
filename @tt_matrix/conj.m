function [b]=conj(a)
%[b]=conj(a)
%complex conjugate of a tt_matrix
b=a;
b.tt=conj(b.tt);