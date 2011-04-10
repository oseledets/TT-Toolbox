function [sz] = size(tt,dim)
%[SZ]=SIZE(TT,DIM)
%Returns the size of the TT-tensor
n=tt.n; m=tt.m;
n0=[n,m];
if exist('dim','var')
    sz = n0(dim);
else
    sz=n0;
    
end
