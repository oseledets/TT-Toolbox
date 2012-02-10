function [y]=col(x, k, varargin)
%Computes the columns of a block TT-tensor 
%   [Y]=COL(X,K)
%   It is useful for block algorithms when r(d+1) is not equal to 1
%   For analogous thing on r(1) see ROW(X,K)
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

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
