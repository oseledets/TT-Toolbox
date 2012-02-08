function [tt] = tucker_to_tt(fc, core)
%Build TT-tensor from the QTT-Tucker representation
%[TT] = TUCKER_TO_TT(FC, CORE) Given QTT-Tucker representation, compute the
%TT-tensor representation 
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

d = core.d;
n = zeros(d,1); rtuck = zeros(d,1);
for i=1:d
    n(i) = size(fc{i},1);
    rtuck(i) = size(fc{i},2);
end;

rcr = core.r;
tt = core;
% keyboard;

for i=1:d
    cr1 = core{i};
    cr1 = permute(cr1, [2, 1,3]);
    cr1 = reshape(cr1, [rtuck(i), rcr(i)*rcr(i+1)]);
    cr1 = fc{i}*cr1; % size n(i), rc1*rc2
    cr1 = reshape(cr1, n(i), rcr(i), rcr(i+1));
    cr1 = permute(cr1, [2, 1, 3]);
    tt{i} = cr1;
end;

end