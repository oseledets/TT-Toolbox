function [tt]=tt_qreshape(tt,s,sz)

% returns .
%
% April 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

d=size(tt,1);
for k=1:d
	szk=size(tt{k});
	szk=[szk,ones(1,s)];
	szkr=szk(s+1:numel(szk));
	tt{k}=reshape(tt{k},[sz(k,:),szkr]);
end

return
end