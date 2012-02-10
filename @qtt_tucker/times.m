function [c]=times(a,b)
%C=A.*B
%   [C]=TIMES(A,B) Hadamard product of two QTT-tuckers
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

if (nargin == 2 )
    % times of tucker factors (whatever it is)
    d = a.dphys;
    c = qtt_tucker;
    c.dphys = d;
    c.tuck = cell(d,1);
    for i=1:d
        c.tuck{i} = a.tuck{i}.*b.tuck{i};
    end;
    c.core = tt_tensor;
    c.core.d = d;
    rca = a.core.r;
    rcb = b.core.r;
    rta = a.core.n;
    rtb = b.core.n;
    c.core.r = rca.*rcb;
    c.core.n = rta.*rtb;
    c.core.ps = cumsum([1; c.core.r(1:d).*c.core.n.*c.core.r(2:d+1)]);
    c.core.core = zeros(c.core.ps(d+1)-1,1);
    for i=1:d
        cra = a.core{i};
        cra = reshape(cra, rca(i)*rta(i)*rca(i+1), 1);
        crb = b.core{i};
        crb = reshape(crb, 1, rcb(i)*rtb(i)*rcb(i+1));
        crc = cra*crb;
        crc = reshape(crc, rca(i), rta(i), rca(i+1), rcb(i), rtb(i), rcb(i+1));
        crc = permute(crc, [1, 4, 2, 5, 3, 6]);
        c.core.core(c.core.ps(i):c.core.ps(i+1)-1) = crc(:);
    end;
end
return
end