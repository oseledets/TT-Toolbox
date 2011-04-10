function [c] = rdivide(a,b)
%[C]=RDIVIDE(A,B)

if isa(a,'tt_tensor') && numel(b) == 1
c=a*(1.0./b);
else
    error('Use rdivide(full(A),full(B)).');
end
