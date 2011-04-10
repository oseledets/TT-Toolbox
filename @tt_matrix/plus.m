function [a]=plus(b,c)
%[A]=PLUS(B,C)
%Adds two TT-matrices B and C
a=tt_matrix; 
a.n=b.n;
a.m=b.m;
a.tt=b.tt+c.tt;
return
end