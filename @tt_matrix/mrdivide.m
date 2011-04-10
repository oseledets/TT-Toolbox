function [c] = mrdivide(a,b)
%[C]=MRDIVIDE(A,B)
%Division of a TT-Matrix by number 
if isa(a,'tt_matrix') && numel(b) == 1
c=a*(1.0./b);
else
    error('Use mrdivide(full(A),full(B)).');
end
