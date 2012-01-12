function [qt]=vector(qt)
% function [qt]=vector(qt)
% Makes a tt_tensor from tt_matrix in factors

d = qt.dphys;
for i=1:d
    qt.tuck{i} = tt_tensor(qt.tuck{i});
end;

end