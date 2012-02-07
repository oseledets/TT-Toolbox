function [reff]=tt_effrank(tt,s)
%Effective rank of the TT-tensor (correct memory definition)
%   [REFF]=TT_EFFRANK(TT,S)
% May 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia

d=size(tt,1);
sz=tt_qsize(tt,s);
n=prod(sz,2);
mem=0;
for k=1:d
	mem=mem+numel(tt{k});
end
a=sum(n(2:d-1));
b=n(1)+n(d);
reff=sqrt(mem/a+b^2/a^2/4)-b/a/2;

return
end
