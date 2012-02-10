function [tt]=tt_random(n,d,r)
%Generates a random tensor
%   [TT]=TT_RANDOM(N,D,R) Generates a random tensor with dimensions specified
%   by (N,D), where N can be a number of an array of dimension D, R is a
%   rank or a array of dimension d+1
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
tt=tt_rand(n,d,r);
return
end
