% function [res] = tt_axpy2(log_a, sign_a, x, log_p, sign_p, y, eps, max_rank)
% Stabilized aX+pY combination, where
% log_a = log(abs(a)), sign_a = sign(a), 
% log_p = log(abs(p)), sign_p = sign(p), 
% eps, max_rank - compression parameters

function [res] = tt_axpy2(log_a, sign_a, x, log_p, sign_p, y, eps, max_rank)

if ((nargin<8)||(isempty(max_rank)))
    max_rank=[];
end;

ax = tt_scal2(x, log_a, sign_a);
py = tt_scal2(y, log_p, sign_p);

res = tt_add(ax,py);
res = tt_compr2(res, eps, max_rank);

end