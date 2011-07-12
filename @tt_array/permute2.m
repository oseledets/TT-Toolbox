function [a] = permute2(a, order)
% function [a] = permute2(a, order)
% Permutes the "physical" dimensions (in each block) of the
% tt_array A.

nsa = a.ns;
tta = a.tt;
ra = tta.r;
d = tta.d;
dp = size(nsa,2);
psa = tta.ps;
cra = tta.core;

for i=1:d
    a1 = cra(psa(i):psa(i+1)-1);
    a1 = reshape(a1, [ra(i), nsa(i,:), ra(i+1)]);
    a1 = permute(a1, [1, (order+1), dp+2]);
    cra(psa(i):psa(i+1)-1) = reshape(a1, ra(i)*prod(nsa(i,:),2)*ra(i+1),1);
end;

tta.core = cra;
a.tt = tta;

end