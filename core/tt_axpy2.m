function [res] = tt_axpy2(log_a, sign_a, x, log_p, sign_p, y, eps, max_rank)
%Returns A*X+P*Y in the TT-format (stabilized version)
%   [RES]=TT_AXPY(A,X,P,Y,EPS,MAX_RANK) Stabilized a*X+p*Y where  log_a = 
%   log(abs(a)), sign_a = sign(a),  log_p = log(abs(p)), sign_p = sign(p), 
%   EPS, MAX_RANK - compression parameters

if ((nargin<8)||(isempty(max_rank)))
    max_rank=[];
end;

ax = tt_scal2(x, log_a, sign_a);
py = tt_scal2(y, log_p, sign_p);

res = tt_add(ax,py);
res = tt_compr2(res, eps, max_rank);

end