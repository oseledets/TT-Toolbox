function [flag] = tt_ismatrix(tt)
%[FLAG]=TT_ISMATRIX(TT)
%Tests if the TT is a matrix or a vector
%
%
% TT Toolbox 1.1, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
if ( (size(tt,1) == 2)) %I suppose it is a matrix in this 2D case 
  flag=false;
elseif (size(size(tt{1}),2)==3 || (size(tt{1},3)==1 && size(tt{1},2) ~= 1 && size(tt{2},3) == 1)) %It is surely a matrix
  flag=true;
else
  flag=false;
end
return
