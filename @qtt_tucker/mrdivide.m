function [c] = mrdivide(a,b)
%[C]=MRDIVIDE(A,B)
%Implement multiplication for a tt_tensor and (in future) 
%matrix-by-vector product

if isa(a,'qtt_tucker') && numel(b) == 1
c=a*(1.0./b);
else
    error('Use mrdivide(full(A),full(B)).');
end
