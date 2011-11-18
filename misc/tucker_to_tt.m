function [tt] = tucker_to_tt(fc, core)
% function [tt] = tucker_to_tt(fc, core)
% Builds the TT tensor from the Tucker representation,
% where fc is the cell array of Tucker factors,
% cr is the Tucker core in TT.

d = core.d;
n = zeros(d,1); rtuck = zeros(d,1);
for i=1:d
    n(i) = size(fc{i},1);
    rtuck(i) = size(fc{i},2);
end;

rcr = core.r;
tt = core;
% keyboard;

for i=1:d
    cr1 = core{i};
    cr1 = permute(cr1, [2, 1,3]);
    cr1 = reshape(cr1, [rtuck(i), rcr(i)*rcr(i+1)]);
    cr1 = fc{i}*cr1; % size n(i), rc1*rc2
    cr1 = reshape(cr1, n(i), rcr(i), rcr(i+1));
    cr1 = permute(cr1, [2, 1, 3]);
    tt{i} = cr1;
end;

end