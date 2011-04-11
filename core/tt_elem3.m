% function [res]=tt_elem3(tt, indices)
% Shrinks physical indices in q-th core of tt, specified by indices{q},
% if size(indices)==size(tt). If size(indices)<size(tt), only last
% size(indices) cores are affected.

function [res]=tt_elem3(tt, indices)

d = size(tt,1);

start_i = d-max(size(indices))+1;
res = tt;
for q=start_i:d
    res{q}=tt{q}(indices{q-start_i+1},:,:);
end;

end