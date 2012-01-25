function [y]=col(x, k, varargin)
%[Y]=COL(X,K)
%Computes the subcolumns of a 
%It is useful for block algorithms when r(d+1) is not equal to 1
%For analogous thing on r(1) see ROW(X,K)
d=x.d;
r=x.r;

if nargin > 2
    for j = 1 : nargin-2
        k = [k, varargin{j}];
    end
end

BlockSize = r(d+1);
if max(k) > BlockSize
    error('TT-tensor.col:Index exceeds dimensions in TT block.');
end

CutMatr = eye(BlockSize);
CutMatr = CutMatr(:, k);
y = x*CutMatr;  

return
end
