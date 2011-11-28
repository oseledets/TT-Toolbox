function [qt]=matrix(qt, sz1, sz2)
% function [qt]=matrix(qt, sz1, sz2)
% Makes factors of qtt_tucker to be tt_matrices
% sz1 and sz2 are cells with n,m dimensions for each factor
% If they are absent, the square matrices are assumed

if (nargin<3)
    sz1 = [];
    sz2 = [];
end;

d = qt.dphys;
for i=1:d
    if (~isempty(sz1))||(~isempty(sz2))
        qt.tuck{i} = tt_matrix(qt.tuck{i}, sz1{i}, sz2{i});
    else
        cursz = qt.tuck{i}.n;
        cursz = sqrt(cursz);
        qt.tuck{i} = tt_matrix(qt.tuck{i}, cursz, cursz);
    end;
end;

end