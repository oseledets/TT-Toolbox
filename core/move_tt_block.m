function [tt] = move_tt_block(tt, spos, epos, eps)
%Performs a bubble movement of a block inside a train
%   [TT] = MOVE_TT_BLOCK(TT, SPOS, EPOS, EPOS) Performs the bubble movement
%   of the SPOS-th block to the position EPOS, the intermediate blocks
%   shift to the vacant places, i.e. by (-1) if spos<eps, and +1,
%   otherwise. Requires truncation with the L2-norm accuracy EPSs
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et. al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

if (spos==epos)
    return;
end;

d = tt.d;
n = tt.n;
r = tt.r;

% QR to spos
for i=1:spos-1
    cr = reshape(tt{i}, r(i)*n(i), r(i+1));
    [cr, rv]=qr(cr, 0);
    cr2 = reshape(tt{i+1}, r(i+1), n(i+1)*r(i+2));
    cr2 = rv*cr2;
    r(i+1) = size(cr, 2);
    tt{i} = reshape(cr, r(i), n(i), r(i+1));
    tt{i+1} = reshape(cr2, r(i+1), n(i+1), r(i+2));
end;
for i=d:-1:spos+1
    cr = reshape(tt{i}, r(i), n(i)*r(i+1));
    [cr, rv]=qr(cr.', 0);
    cr2 = reshape(tt{i-1}, r(i-1)*n(i-1), r(i));
    cr2 = cr2*(rv.');
    r(i) = size(cr, 2);
    tt{i} = reshape(cr.', r(i), n(i), r(i+1));
    tt{i-1} = reshape(cr2, r(i-1), n(i-1), r(i));
end;

% Now, start permutation
if (spos<epos) % From left to right  
    for i=spos:epos-1
        cr1 = reshape(tt{i}, r(i)*n(i), r(i+1));
        cr2 = reshape(tt{i+1}, r(i+1), n(i+1)*r(i+2));
        cr = cr1*cr2;
        cr = reshape(cr, r(i), n(i), n(i+1), r(i+2));
        cr = permute(cr, [1, 3, 2, 4]);
        tempvar = n(i); n(i) = n(i+1); n(i+1) = tempvar;
        cr = reshape(cr, r(i)*n(i), n(i+1)*r(i+2));
        [u,s,v]=svd(cr, 'econ');
        s = diag(s);
        nrm = norm(s);
        r(i+1) = my_chop2(s, eps*nrm/sqrt(d-1));
        tt{i} = reshape(u(:,1:r(i+1)), r(i), n(i), r(i+1));
        v = diag(s(1:r(i+1)))*(v(:,1:r(i+1))');
        tt{i+1} = reshape(v, r(i+1), n(i+1), r(i+2));
    end;
end;

if (epos<spos) % From right  to left
    for i=spos:-1:epos+1
        cr1 = reshape(tt{i}, r(i), n(i)*r(i+1));
        cr2 = reshape(tt{i-1}, r(i-1)*n(i-1), r(i));
        cr = cr2*cr1;
        cr = reshape(cr, r(i-1), n(i-1), n(i), r(i+1));
        cr = permute(cr, [1, 3, 2, 4]);
        tempvar = n(i); n(i) = n(i-1); n(i-1) = tempvar;
        cr = reshape(cr, r(i-1)*n(i-1), n(i)*r(i+1));
        [u,s,v]=svd(cr, 'econ');
        s = diag(s);
        nrm = norm(s);
        r(i) = my_chop2(s, eps*nrm/sqrt(d-1));
        tt{i} = reshape(v(:,1:r(i))', r(i), n(i), r(i+1));
        u = u(:,1:r(i))*diag(s(1:r(i)));
        tt{i-1} = reshape(u, r(i-1), n(i-1), r(i));
    end;
end;

end
