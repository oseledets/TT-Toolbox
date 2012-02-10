function [r]=rank(a,varargin)
%All or individual TT-ranks of a TT-matrix
%   [R]=RANK(A) Returns all (d+1) ranks of a TT-matrix A
%   
%   [R]=RANK(A,IND) Return the TT-rank with the index IND
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
if (nargin==1)
  r=a.tt.r;
else
  i=varargin{1};
  r=rank(a.tt,i);
end
return
end