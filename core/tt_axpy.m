function [res]=tt_axpy(a, x, p, y, eps, max_rank)
%Returns A*X+P*Y in the TT-format. 
%   [RES]=TT_AXPY(A,X,P,Y,EPS,MAX_RANK) Returns a*X+p*Y in TT format. 
%   Compress each summand with respect to its norm. Compress the result to 
%   accuracy EPS and (optional) maximal rank MAX_RANK
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
% d = size(x, 1);
if ((nargin<6)||(isempty(max_rank)))
    max_rank=[];
end;

% nrmx = tt_dot(x,x);
% nrmy = tt_dot(y,y);
% eps1 = (abs(a)*nrmx+abs(p)*nrmy)/(abs(a)*nrmx);
% eps2 = (abs(a)*nrmx+abs(p)*nrmy)/(abs(p)*nrmy);
% 
% % aa = abs(a)
% % ap = abs(p)
% 
% if (eps1>2)
% %     eps1
%     x = tt_compr2(x, min(eps*eps1, 0.5), max_rank);
% end;
% if (eps2>2)
% %     eps2
%     y = tt_compr2(y, min(eps*eps2, 0.5), max_rank);
% end;

ax = tt_scal(x, a);
py = tt_scal(y, p);

res = tt_add(ax,py);
res = tt_compr2(res, eps, max_rank);

% res = tt_stabilize(res,0);

end