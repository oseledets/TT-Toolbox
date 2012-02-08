function [fc,core]=qtt_tucker_m(tt, sz, tol)
%Compute the QTT-Tucker representation of TT
%   [FC,CORE]=QTT_TUCKER_M(TT, SZ, TOL) Computes the QTT-Tucker 
%   representation of TT-tensor. This is a technical function: use 
%   corresponding constructors of the QTT-Tucker class
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


d = tt.d;
n = tt.n;
rcr = tt.r;

fc = cell(d,1);

% Initial orthogonalization - 1->d
for i=1:d-1
    cr = tt{i};
    cr = reshape(cr, rcr(i)*n(i), rcr(i+1));
    [cr, rv]=qr(cr, 0);
    cr2 = tt{i+1};
    cr2 = rv*reshape(cr2, rcr(i+1), n(i+1)*rcr(i+2));
    rcr(i+1) = size(cr,2);
    tt{i} = reshape(cr, rcr(i), n(i), rcr(i+1));
    tt{i+1} = reshape(cr2, rcr(i+1), n(i+1), rcr(i+2));
end;

% Split tucker factors from d->1
core = tt;
nextcr = core{d};
rtuck = zeros(d,1);
for i=d:-1:1
    cr = nextcr;
    cr = permute(cr, [2, 1, 3]);
    cr = reshape(cr, n(i), rcr(i)*rcr(i+1));
    [u,s,v]=svd(cr, 'econ');
    s = diag(s);
    nrm = norm(s);
    r = my_chop2(s, tol*nrm); % !!!!!!
    rtuck(i) = r;
    if (i>1)
        fc{i} = u(:,1:r);
        cr = diag(s(1:r))*v(:,1:r)';
        cr = reshape(cr, r, rcr(i), rcr(i+1));
        cr = permute(cr, [1, 3, 2]);
        cr = reshape(cr, r*rcr(i+1), rcr(i));
        [cr,rv] = qr(cr, 0);
        cr2 = core{i-1};
        cr2 = reshape(cr2, rcr(i-1)*n(i-1), rcr(i));
        cr2 = cr2*(rv.');
        rcr(i) = size(cr, 2);
        nextcr = reshape(cr2, rcr(i-1), n(i-1), rcr(i));
        cr = reshape(cr, r, rcr(i+1), rcr(i));
        core{i} = permute(cr, [3, 1, 2]);
    else
        cr = v(:,1:r)';
        cr = reshape(cr, r, rcr(i), rcr(i+1));
        core{i} = permute(cr, [2, 1, 3]);
        fc{i} = u(:,1:r)*diag(s(1:r));
    end;    
end;

% Quantics approximation
for i=1:d
    %d0 = log2(n(i));
    fc{i} = tt_tensor(fc{i});
    szc=sz{i}; d0=numel(szc); m0=szc(end); szc=szc(:);
    %fc{i} = tt_reshape(fc{i}, [2*ones(1,d0-1), 2*rtuck(i)], tol);
    fc{i}=tt_reshape(fc{i},[szc(1:end-1); m0*rtuck(i)],tol);
    nocore = fc{i}{d0}; % rqtt, 2*rtuck, 1
    rqtt = size(nocore, 1);
    if (i<d)
        % Play with QRs: fc{i}->cr{i}->cr{i+1}->fc{i+1}        
        nocore = reshape(nocore, rqtt*m0, rtuck(i));
        [ocore, nocore]=qr(nocore, 0);               
        cr1 = core{i};
        cr1 = permute(cr1, [2, 1, 3]);
        cr1 = reshape(cr1, rtuck(i), rcr(i)*rcr(i+1));
        cr1 = nocore*cr1;
        rtuck(i) = size(cr1, 1);
        fc{i}{d0} = reshape(ocore, rqtt, m0, rtuck(i));
        
        cr1 = reshape(cr1, rtuck(i)*rcr(i), rcr(i+1));
        [cr1, rv]=qr(cr1, 0);
        cr2 = core{i+1};
        cr2 = reshape(cr2, rcr(i+1), rtuck(i+1)*rcr(i+2));
        cr2 = rv*cr2;
        rcr(i+1) = size(cr1, 2);
        cr1 = reshape(cr1, rtuck(i), rcr(i), rcr(i+1));
        core{i} = permute(cr1, [2,1,3]);
        cr2 = reshape(cr2, rcr(i+1), rtuck(i+1), rcr(i+2));
        cr2 = permute(cr2, [1, 3, 2]);
        cr2 = reshape(cr2, rcr(i+1)*rcr(i+2), rtuck(i+1));
        [cr2, rv]=qr(cr2, 0);
        fc{i+1} = fc{i+1}*(rv.');
        rtuck(i+1) = size(cr2, 2);
        cr2 = reshape(cr2, rcr(i+1), rcr(i+2), rtuck(i+1));
        core{i+1} = permute(cr2, [1, 3, 2]); 
    else
        fc{i}{d0} = reshape(nocore, rqtt, m0, rtuck(i));
    end;
end;

end