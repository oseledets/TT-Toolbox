% function [res] = full_to_qtt(vec, eps, d)
% Compress full tensor vec to QTT tensor with accuracy eps. If d is specified,
% set the amount of data to 2^d. Otherwise, vec is treated as vector of
% size n-by-1

function [res] = full_to_qtt(vec, eps, d)

N = max(size(vec));
if ((nargin<3)||(isempty(d)))
d = round(log2(N));
end;

sz = 2*ones(1, d);

res = reshape(vec, sz);

res = full_to_tt(res, eps);
end