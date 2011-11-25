function [c] = mtimes(a,b)
%[C]=MTIMES(A,B)
%Elementwise multiplication for the QTT_TUCKER class

if isa(a,'double') && isa(b,'qtt_tucker')
    c=b;
    c.core=a*b.core;
elseif isa(a,'qtt_tucker') && isa(b,'double') 
    c=a;
    c.core=a.core*b;
    
else
    error('Use mtimes(full(A),full(B)).');
end
