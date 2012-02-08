function [tt_res] = tt_add(tt1,tt2)
%Adds two TT-tensors in TT1 format.
%   [TT_RES]=TT_ADD(TT1,TT2) - adds two tensors in TT1.0 format. Please
%   avoid the usage: it will be removed from future releases. Use the
%   object-oriented version.
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



d=size(tt1,1);
tt_res=cell(d,1);
tt_res{1}=[tt1{1} tt2{1}];
tt_res{d}=[tt1{d} tt2{d}];

for i=2:d-1
    core1=tt1{i};
    core2=tt2{i};
    ncur=size(core1,1);
    r21=size(core1,2);
    r31=size(core1,3);
    r22=size(core2,2);
    r32=size(core2,3);
    core=zeros(ncur,r21+r22,r31+r32);
    core(1:ncur,1:r21,1:r31)=core1;
    core(1:ncur,(r21+1):(r21+r22),(r31+1):(r31+r32))=core2;
    tt_res{i}=core;
end
return
end
