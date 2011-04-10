% function [res]=qtt_to_full(qtt,ns)
% Deploys QTT tensor in full format, specified by sizes in ns.
% default sizes is 2^d-by-1

function [res]=qtt_to_full(qtt,ns)

d = size(qtt,1);
res = tt_to_full(qtt);

if ((nargin<2)||(isempty(ns))) ns=[2^d, 1]; end;

res = reshape(res, ns);

end