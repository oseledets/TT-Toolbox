function [nrm]=normest(t)
%2nd matrix norm of the TT-matrix
%   [nrm]=normest(t)
%   Fast 2-norm estimation by 5 power iterations.
%   The matrix should be square
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

v = tt_rand(t.n, t.tt.d, 5);
for i=1:5
    v = v/norm(v);
    w = t*v;
    w = round(w, 1e-3);
    nrm = norm(w);
    v = w;
end;
end