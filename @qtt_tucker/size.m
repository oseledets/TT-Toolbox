function [sz] = size(tt,dim)
%[SZ]=SIZE(TT,DIM)
%Returns the size of the TT-tensor

if exist('dim','var')
    sz = tt.n(dim);
else
    sz=tt.n;
    
end
