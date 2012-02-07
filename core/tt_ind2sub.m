function [ind] = tt_ind2sub(siz,ndx)
%Correct conversion of an index to a multiindex
%   [IND]=TT_IND2SUB(SIZ,NDX) 
%
%
% TT Toolbox 2.1, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
%ind=tt_ind2sub(size,ndx)

%nout = max(nargout,1);
nout = size(siz,2);

siz = double(siz);
ind=zeros(1,nout);
if length(siz)<=nout,
  siz = [siz ones(1,nout-length(siz))];
else
  siz = [siz(1:nout-1) prod(siz(nout:end))];
end
n = length(siz);
k = [1 cumprod(siz(1:end-1))];
for i = n:-1:1,
  vi = rem(ndx-1, k(i)) + 1;         
  vj = (ndx - vi)/k(i) + 1; 
  ind(i) = vj; 
  ndx = vi;     
end
