function [M]=IpaS(d, a)
%A special bidiagonal matrix in the QTT-format
%   [M]=IPAS(D, A)
%   Generates I+a*S_{-1} matrix in the QTT-format:
%   1 0 0 0
%   a 1 0 0
%   0 a 1 0
%   0 0 a 1
% Convenient for Crank-Nicolson and time gradient matrices
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


M = cell(d,1);
M{1} = zeros(2,2,0);
M{1}(:,:,1)=[1,0;a,1];
if (d>1)
    M{1}(:,:,2)=[0,a;0,0];
end;
for i=2:d-1
    M{i}=zeros(2,2,2,2);
    M{i}(:,:,1,1)=eye(2);
    M{i}(:,:,2,1)=[0,0;1,0];
    M{i}(:,:,2,2)=[0,1;0,0];
end;
if (d>1)
    M{d} = zeros(2,2,2);
    M{d}(:,:,1)=eye(2);
    M{d}(:,:,2)=[0,0;1,0];
end;
M=tt_matrix(M); %Bydlocode
end
