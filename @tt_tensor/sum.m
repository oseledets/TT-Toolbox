function [val] = sum(tt,chunk_start,chunk_end)
%[V] = SUM(TT)
%[V] = SUM(TT, dim) 
%[V] = SUM(TT, chunk_start,chunk_end) 

%[V] = SUM(TT)
%   Computes the sum of all entries of a TT-tensor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[V] = SUM(TT, dim) 
% dim: integer,  1 <= dim <= tt.d
%   Computes the sums of entries along the dimension indicated by "dim",
%   similar to Matlab built-in function SUM(tensor, dim);
% return a TT tensor of order (tt.d-1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[V] = SUM(TT, chunk_start,chunk_end) 

%chunk_start: integer,  1 <= chunk_start <= tt.d
%chunk_end: integer, chunk_start <= chunk_end <= tt.d

%   Computes the sums of entries of the sub-tensors indicated by
%   chunk_start and chunk_end. The sub-tensor taken is the same as in
%   chunk(tt, chunk_start,chunk_end).
%   return a TT tensor of order (tt.d- (chunk_end -chunk_start +1) ).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    e = tt_ones(size(tt) );
    val = dot(e,tt);
else
    if (nargin==2)
     chunk_end  = chunk_start;   
    end
    e = tt_ones(tt.n(chunk_start:chunk_end));
    val = dot(e,tt, chunk_start,chunk_end);   
end 
end
