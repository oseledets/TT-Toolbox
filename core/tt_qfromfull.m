function [tt]=tt_qfromfull(full,s,d,eps)

% Approximates a full-format s-dimensional 2^{d} x...x 2^{d}-tensor
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


full=reshape(full,2*ones(1,s*d));

prm=zeros(1,s*d);
for k=1:d
	for l=1:s
		prm(s*k+1-l)=k+(s-l)*d;
	end
end

full=permute(full,prm);
full=reshape(full,2^s*ones(1,d));

tt=full_to_tt(full,eps);

for k=1:d
	[n,p,q]=size(tt{k});
	tt{k}=reshape(tt{k},[2*ones(1,s),p,q]);
end

return
end