function y=horzcat(varargin)
%% TT-Toolbox 2.2.1, 2009-2016
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel

%for questions on this particular code email Alexey Boyko, Skolkovo Institute of Science and Technology,
% a.boyko@skoltech.ru or alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2016
%%

%now you can assemble rectangular matrices horizontally like 
%ttAA=[ttA ttA]

%Attention!
%Due to limitations of MATLAB, for simultaneous horizontal and vertical concatenation 
% use
% ttM=tt_blockwise({ttA ttB ttC; ttD ttE ttF ;ttG ttH ttI }) 
% instead of 
% ttM=[ttA ttB ttC; ttD ttE ttF ;ttG ttH ttI ]

y=tt_blockwise(varargin);
end