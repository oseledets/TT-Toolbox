function [qt]=matrix(qt, sz1, sz2)
%Converts a vector to a matrix in the QTT-Tucker format
%   [QT]=MATRIX(QT, SZ1, SZ2) Makes factors of QTT-Tucker to be TT-matrices
%   SZ1 and SZ2 are cell arrays with N, M dimensions for each factor
%   If they are absent, the square matrices are assumed
% 
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

if (nargin<3)
    sz1 = [];
    sz2 = [];
end;

d = qt.dphys;
for i=1:d
    if (~isempty(sz1))||(~isempty(sz2))
        qt.tuck{i} = tt_matrix(qt.tuck{i}, sz1{i}, sz2{i});
    else
        cursz = qt.tuck{i}.n;
        cursz = sqrt(cursz);
        qt.tuck{i} = tt_matrix(qt.tuck{i}, cursz, cursz);
    end;
end;

end