function [ful]=tt_qtofull(tt,s)

% Returns a tensor _full_, given in a QTT decomposition _tt_,
% k-th core of which is indexed by
% i^{k}_{1} ,..., i^{k}_{s},
% k=1,...,d,
% in the full-format representation.
% The output tensor _full_ is indexed by
% i^{1}_{1} ,..., i^{d}_{1} ,...,..., i^{1}_{s} ,..., i^{d}_{s},
%
% April 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

d=size(tt,1);
sz=tt_qsize(tt,s);

tt=tt_qreshape(tt,s,prod(sz,2));
ful=full(tt_tensor(tt));

ful=reshape(ful,reshape(sz',1,numel(sz)));

prm=zeros(1,s*d);
for k=1:d
	for l=1:s
		prm(s*k+1-l)=k+(s-l)*d;
	end
end

ful=ipermute(ful,prm);
ful=reshape(ful,[prod(sz,1),1]);

return
end