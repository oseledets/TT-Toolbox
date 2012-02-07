function [tt]=tt_kron7(tt1,tt2,s)
%Kronecker product of two TT-tensors in the reverse order
%   [TT]=TT_KRON7(TT1,TT2,S)
% May 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia

d1 = size(tt1,1);
d2 = size(tt2,1);

tt1{1}=permute(tt1{1},[(1:s),s+2,s+1]);
%tt2{1}=permute(tt2{1},[(1:s),s+2,s+1]);

tt = cell(d1+d2,1);
for i=1:d2
	tt{i}=tt2{i};
end
for i=1:d1
	tt{d2+i}=tt1{i};
end

%tt{1}=permute(tt{1},[(1:s),s+2,s+1]);

return
end
