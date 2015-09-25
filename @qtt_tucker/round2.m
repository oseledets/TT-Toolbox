function [tt]=round2(tt,varargin)
%Approximate QTT-Tucker with another one with specified accuracy, ver 2.
%   [QTT]=ROUND2(QTT,EPS) Approximate QTT-Tucker tensor with relative 
%   accuracy EPS
%
%   [QTT]=ROUND2(QTT,EPS,RMAX) Approximate QTT-Tucker tensor with relative
%   accuracy 
%   EPS and maximal rank RMAX. RMAX can be array of ranks or a number
%
% Thanks to Simon Etter (ETH Zurich) for correcting a mistake in the
% original round algorithm that led to a factor d in the truncation error
% instead of sqrt(d) in this new version.
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com

d=tt.dphys;
core=tt.core;
tuck=tt.tuck;
eps=varargin{1};
tolcorr = 0; % number of all SVDs, to adjust eps by sqrt(tolcorr)
for i=1:d
    tolcorr = tolcorr+tuck{i}.d;
end;
tolcorr = tolcorr + d-1;
rmax = [];
if (nargin==3)
    rmax = varargin{2};
end;
% If we have a qtt_tucker matrix, remember its sizes and stretch into the
% vector
ismatrix = 0;
if (isa(tuck{1}, 'tt_matrix'))
    ismatrix = 1;
    curn = cell(d,1);
    curm = cell(d,1);
end;
for i=2:d
    if (ismatrix)
        curn{i} = tuck{i}.n;
        curm{i} = tuck{i}.m;
        tuck{i} = tt_tensor(tuck{i});
    end;
   [tuck{i},rm]=qr(tuck{i},'lr');
   core{i}=ten_conv(core{i},2,rm.');
end
% All factors except the first one are left-orthogonal, orthogonalize the core now
[core,rm]=qr(core, 'rl');
core{1} = rm*core{1};

% Sweep over physical dimensions
for i=1:d
    % Orthogonalize RL and truncate LR factors
    L = tuck{i}.d;
    tucki = core2cell(tuck{i});
    tucki = [tucki; {core{i}}];
    [rc1,rt,rc2]=size(tucki{L+1});
    tucki{L+1} = permute(tucki{L+1}, [2,1,3]);
    tucki{L+1} = reshape(tucki{L+1}, rt, rc1*rc2);
    r1 = rt;
    n1 = rc1*rc2;
    r2 = 1;
    % Orth
    for j=L+1:-1:2
        [u,rm]=qr(tucki{j}.', 0);
        rnew = size(u,2);
        [r0,n0,~]=size(tucki{j-1});
        tucki{j-1} = reshape(tucki{j-1}, r0*n0, r1);
        tucki{j-1} = tucki{j-1}*rm.';
        tucki{j-1} = reshape(tucki{j-1}, r0, n0*rnew);
        tucki{j} = reshape(u.', rnew, n1, r2);
        r1 = r0;
        n1 = n0;
        r2 = rnew;
    end;
    % Trunc
    for j=1:L
        tucki{j} = reshape(tucki{j}, r1*n1, r2);
        [u,s,v]=svd(tucki{j}, 'econ');
        s = diag(s);
        rnew = my_chop2(s, eps*norm(s)/sqrt(tolcorr));
        if (~isempty(rmax))
            rnew = min(rnew, rmax);
        end;        
        u = u(:,1:rnew);
        v = diag(s(1:rnew))*v(:,1:rnew)';
        tucki{j} = reshape(u, r1, n1, rnew);
        [r2,n2,r3]=size(tucki{j+1});
        tucki{j+1} = reshape(tucki{j+1}, r2, n2*r3);
        tucki{j+1} = v*tucki{j+1};
        tucki{j+1} = reshape(tucki{j+1}, rnew, n2, r3);
        r1 = rnew;
        n1 = n2;
        r2 = r3;
    end;
    rt = size(tucki{L+1}, 1);
    tuck{i} = cell2core(tt_tensor, tucki(1:L));
    % Truncate the core block
    tucki{L+1} = reshape(tucki{L+1}, rt, rc1, rc2);
    tucki{L+1} = permute(tucki{L+1}, [2,1,3]);
    if (i<d)
        tucki{L+1} = reshape(tucki{L+1}, rc1*rt, rc2);
        [u,s,v]=svd(tucki{L+1}, 'econ');
        s = diag(s);
        rnew = my_chop2(s, eps*norm(s)/sqrt(tolcorr));
        if (~isempty(rmax))
            rnew = min(rnew, rmax);
        end;
        tucki{L+1} = u(:, 1:rnew);
        tucki{L+1} = reshape(tucki{L+1}, rc1, rt, rnew);
        v = diag(s(1:rnew))*v(:,1:rnew)';
        core2 = core{i+1};
        [~,rt2,rc3]=size(core2);
        core2 = reshape(core2, rc2, rt2*rc3);
        core2 = v*core2;
        core{i+1} = reshape(core2, rnew, rt2, rc3);
    end;
    core{i} = tucki{L+1};
end;

tt.core = core;
tt.tuck = tuck;
end