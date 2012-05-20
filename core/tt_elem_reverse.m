function [tt]=tt_elem_reverse(tt)
% Elementwise reverse of a tensor train vector
% NOT tested on a multidimensional case
%
% use tt_qutrtoepl(core(tt_shf(d)'*tt_elem_reverse(tt))) - to build
% upper-toeplitz matrix in "normal", blin, indexing, i.e. second element
% goes to S^{+1}, usw.

d = tt.d;
tt = core2cell(tt);

for i=1:d
    n = size(tt{i}, 2);
    tt{i} = tt{i}(:,(n:-1:1),:);
end;

tt = cell2core(tt_tensor, tt);

end