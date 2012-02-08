function [mem]=tt_mem(tt)
% Memory required to store the tensor TT
%   [MEM]=TT_MEM(TT) Computes the total number of elements in the TT-tensor
%   in the TT1.0 format. Please avoid its usage: it will be removed in
%   future releases. Use mem() from the object-oriented version.
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
mem=0;
d=size(tt,1);
for k=1:d
    mem=mem+numel(tt{k});
end
end
