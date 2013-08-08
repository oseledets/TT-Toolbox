function [b]=real(a)
%Compute a real part of TT-matrix 
%   [B]=REAL(A)
%
% TT-Toolbox 2.2, 2009-2013
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%
% For details see Thm4 in Dolgov, Khoromskij, Savostyanov,
% "Superfast Fourier transform using the QTT format"
% http://dx.doi.org/10.1007/s00041-012-9227-4
%
%---------------------------

b=a;
b.tt=real(b.tt);
return
end