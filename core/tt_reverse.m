function [tt2]=tt_reverse(tt,phd)
% Returns a reversed TT representation
%   [TT2]=TT_REVERSE(TT,PHD)
% January 25, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

d=size(tt,1);
tt2=cell(d,1);

tt2{1}=tt{d};
for k=2:d-1
	prm=(1:phd+2);
	prm(phd+1)=phd+2;
	prm(phd+2)=phd+1;
	tt2{k}=permute(tt{d+1-k},prm);
end
tt2{d}=tt{1};

return
end