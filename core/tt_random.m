function [tt]=tt_random(n,d,r)
%[TT]=TT_RANDOM(N,D,R)
%Generates a random D-dimensional tensor TT with mode sizes N and ranks R
%N can be either mode vector or a number
%R can be either mode vector or a number
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
tt=tt_rand(n,d,r);
return
end
