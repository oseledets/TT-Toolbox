function [tt]=qtttucker_to_tt(fc, cr)
% function [tt]=qtttucker_to_tt(fc, cr)
% Builds the TT tensor from the QTT-Tucker representation,
% where fc is the Tucker factors in QTT (cell array of @tt_tensor),
% cr is the Tucker core in TT.

d = cr.d;
for i=1:d
    fc{i} = tt_reshape(fc{i}, prod(fc{i}.n), [], 1, fc{i}.r(fc{i}.d+1));
    fc{i} = fc{i}{1};
    fc{i} = reshape(fc{i}, size(fc{i}, 2), size(fc{i}, 3));
end;
tt = tucker_to_tt(fc, cr);

end