function [mem] = tt_save(filename, tt, precision)
%Saves TT-tensors into a Fortran-readable SDV format
%   [MEM] = TT_SAVE(FILENAME, TT, [PRECISION]) Saves a TT-tensor  to a 
%   filename in SDV =) format, precision: 0 - double, 1 - single  returns 
%   the bytes written to a file, -1 in case of errors.
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

if (nargin<3)||(isempty(precision))
    precision = 0; % double
end;

mem = tt_write(filename, tt.d, tt.r, tt.n, tt.core, precision);

end
