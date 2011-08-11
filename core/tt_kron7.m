function [tt]=tt_kron7(tt1,tt2)

% May 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia

d1 = size(tt1,1);
d2 = size(tt2,1);
tt = cell(d1+d2,1);
for i=1:d2
	tt{i}=tt2{i};
end
tt{d2+1}=permute(tt1{1}, [1,3,2]);
for i=2:d1
	tt{d2+i}=tt1{i};
end

return
end
