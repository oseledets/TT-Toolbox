function [H]=household(v)

vm = tt_matrix(v,2,1);
nrm = dot(v,v);
H = tt_matrix(tt_eye(2,v.d)) - (vm*vm')*(2/nrm);

end