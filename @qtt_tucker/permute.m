function [tt1] = permute(tt,prm,varargin)
%[TT1]=PERMUTE(TT,PRM)
%[TT1]=PERMUTE(TT,PRM,EPS)
%PERMUTE Permute dimensions for a ttensor.
%It is a difficult procedure in the TT-format, and actually would
%require approximate computation, thus auxiliary parameter EPS
%can be specified. You will be warned, if the permution is extremely
%rank-increasing
if ( nargin == 2 )
elseif ( nargin == 3 )
end
error('Permute function is not yet implemented');
return
end



