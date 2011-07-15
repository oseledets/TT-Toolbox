function [a] = full(tt, nfull)
%[A]=FULL(TT,[nfull])
%Creates full tensor from the array
%Beware of OOM for the full representation!

ns = tt.ns;
d0 = size(ns,1);
dp = size(ns,2);
if (nargin<2)||(isempty(nfull))
    nfull = prod(ns, 1);
end;
if (numel(nfull)==1)
    nfull=[nfull,1];
end;
tt1=tt.tt;
a=full(tt1);

allsizes = reshape(ns', 1, dp*d0);
a = reshape(a, allsizes); % now it is i1,j1,k1,i2,j2,k2,...

perm = ones(1,dp*d0);

for j=1:dp
    for i=1:d0
        perm(i+(j-1)*d0) = j+(i-1)*dp;
    end;
end;

a = permute(a, perm);
a = reshape(a, nfull);

return;
