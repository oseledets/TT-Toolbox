function [r]=rank(a,varargin)
%[R]=RANK(A)
%[R]=RANK(A,IND)
%Returns either all ranks of the TT-tensor A, or rank with index IND
if (nargin==1)
  r=a.r;
else
  i=varargin{1};
   r=a.r(i);
end
return
end