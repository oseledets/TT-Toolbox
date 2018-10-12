function [ind]=tt_sub2ind(siz, idx)
%Correct conversion of a multiindex to an index
%   [IND]=TT_SUB2IND(SIZ,IDX) 
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

% Make it vectorized
n=numel(siz);
k = reshape(siz, n, 1);
k(1)=1;
k(2:n)=siz(1:n-1);
k = cumprod(k);
ind = idx-1;

% Apply to blocks of the second dimension (useful for QTT)
I = size(ind,1);
ind = reshape(ind, I*n, []);
ind = ind.';
ind = reshape(ind, [], n);

% Merging of multi-indices
ind = ind*k;

ind = reshape(ind, [], I);
ind = ind.';

ind = ind+1;

end
