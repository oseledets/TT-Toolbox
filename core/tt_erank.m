function [rk]=tt_erank(tt)
%[RK]=TT_ERANK(TT)
%Computes the effective rank of the TT tensor,
%assuming all mode sizes are equal
%Result is RK
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
  rk = sqrt(tt_mem(tt)/(sum(tt_size(tt))));
 
