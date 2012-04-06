function [tt] = qtttucker_to_linqtt(qt, eps)
%Compute the standard QTT-representation from the QTT-Tucker
%   [TT]=QTTTUCKER_TO_LINQTT(TT, EPS) Computes the "linear" QTT
%   representation from the QTT-TUCKER. This is a technical function: use
%   the corresponding constructor of the QTT-Tucker class.
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


dphys = qt.dphys;
tuck = qt.tuck;
core = qt.core;
rc = core.r;
rt = core.n;
ns = cell(dphys,1);
for i=1:dphys
    ns{i} = tuck{i}.n;    
end;

tt = []; % we will kron the chunks
%  tt = tt_tensor;
%  tt.d = sum(L);
%  tt.ps = ones(tt.d+1,1);
%  tt.n = zeros(tt.d,1);
%  tt.r = zeros(tt.d+1,1);

% Move the nonorth to the first factor
for i=1:dphys
    [tuck{i},rv]=qr(tuck{i}, 'lr');
    cr = core{i};
    cr = permute(cr, [2,1,3]);
    cr = reshape(cr, rt(i), rc(i)*rc(i+1));
    cr = rv*cr;
    rt(i) = size(rv, 1);
    cr = reshape(cr, rt(i), rc(i), rc(i+1));
    core{i} = permute(cr, [2, 1, 3]);
end;
[core,nrm]=qr(core, 'rl');
cr = core{1};
cr = nrm*cr;
core{1} = cr;
rc = core.r;
%  cr = permute(cr, [1,3,2]);
%  cr = reshape(cr, rc(1)*rc(2), rt(1));
%  [q,rv]=qr(cr, 0);
%  tuck{1} = tuck{1}*(rv.');
%  rt(1) = size(q,2);
%  q = reshape(q, rc(1), rc(2), rt(1));
%  core{1} = permute(q, [1, 3, 2]);

% Now, we can start to move the ranks to the linear model
for i=1:dphys
    curtuck = tuck{i};
    cr = core{i};
    cr = permute(cr, [2, 1, 3]);
    cr = reshape(cr, rt(i), rc(i)*rc(i+1));
    curtuck = curtuck*cr; % now the last rank is rc(i)*rc(i+1)
    curtuck = tt_reshape(curtuck, [ns{i}; rc(i); rc(i+1)], eps);
    % Now the most dangerous: the movement of the block rc(i) to the beginning
    curtuck = move_tt_block(curtuck, numel(ns{i})+1, 1, eps);
    if (numel(ns{i})>1)
        curtuck = tt_reshape(curtuck, [rc(i)*ns{i}(1); ns{i}(2:end-1); ns{i}(end)*rc(i+1)], eps);
    else
        curtuck = tt_reshape(curtuck, rc(i)*ns{i}(1)*rc(i+1), eps);
    end;
    % the data is in fact done; update the rank and n formally
    curtuck.r(1) = rc(i);
    curtuck.n(1) = ns{i}(1);
    curtuck.r(numel(ns{i})+1) = rc(i+1);
    curtuck.n(numel(ns{i})) = ns{i}(numel(ns{i}));    
    % Ha-ha, that's not all; QR to the rest of the tucker core, so we will be able to proceed further
    if (i<dphys)
	[curtuck, rv]=qr(curtuck, 'lr');
	cr2 = core{i+1};
	cr2 = reshape(cr2, rc(i+1), rt(i+1)*rc(i+2));
	cr2 = rv*cr2;
	rc(i+1) = size(rv, 1);
	core{i+1} = reshape(cr2, rc(i+1), rt(i+1), rc(i+2));
    end;
    % Now, stuff to the TT
    tt = kron(tt, curtuck);
end;


end
