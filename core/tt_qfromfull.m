function [tt]=tt_qfromfull(ff,s,d,eps,q)

% Approximates a full-format s-dimensional q^{d} x...x q^{d}-tensor
% full, indexed by 
% i^{1}_{1} ,..., i^{d}_{1} ,...,..., i^{1}_{s} ,..., i^{d}_{s},
% by a QTT decomposition _tt_,
% k-th core of which is indexed by
% i^{k}_{1} ,..., i^{k}_{s},
% k=1,...,d
%
% August 11, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

if (nargin<4) || isempty (eps)
    eps=1.e-8;
end

if (nargin<5) || isempty (q)
    q=2;
end

ff=reshape(ff,[q*ones(1,s*d),1]);

prm = 1:s*d; prm = reshape(prm,[d,s]); prm = prm'; prm = reshape(prm,[1,d*s]);

ff=permute(ff,[prm,1+s*d]);
ff=reshape(ff,[q^s*ones(1,d),1]);

tt=core(tt_tensor(ff,eps));
if (numel(tt) > d)
	tt{d}=tt{d}*tt{d+1};
	tt=tt(1:d);
end

for k=1:d
	[~,p,r]=size(tt{k});
	tt{k}=reshape(tt{k},[q*ones(1,s),p,r]);
end

return
end
