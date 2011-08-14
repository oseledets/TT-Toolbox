function [sz]=tt_qsize(tt,s)

% For a QTT decomposition _tt_,
% k-th core of which is indexed by
% i^{k}_{1} ,..., i^{k}_{s},
% k=1,...,d,
% returns the array of mode lengths
% sz = [n^{1}_{1} ,..., n^{1}_{s};
%                 ,...,
%       n^{d}_{1} ,..., n^{d}_{s}]
%
% April 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

d=size(tt,1);
sz=zeros(d,s);

for k=1:d
	szk=size(tt{k});
	szk=[szk,ones(1,s)];
	sz(k,1:s)=szk(1:s);
end

return
end