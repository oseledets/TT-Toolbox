function [tt]=tt_shf(d)
%Upper shift matrix in the QTT-format
%   [TT]=TT_SHF(D) computes an upper shift matrix of size 2^D x 2^D in the
%   TT-format
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et. al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

tt=cell(d,1);
mat1=eye(2); mat2=zeros(2); mat3=zeros(2);
mat2(1,2)=1; mat3(2,1)=1;

tt{1}=zeros(2,2,2);
tt{1}(:,:,1)=mat2; tt{1}(:,:,2)=mat3;
tt{d}=zeros(2,2,2);
tt{d}(:,:,1)=mat1; 
tt{d}(:,:,2)=mat2;
for k=2:(d-1)
   tt{k}=zeros(2,2,2,2);
   tt{k}(:,:,1,1)=mat1;
   tt{k}(:,:,2,1)=mat2;
   tt{k}(:,:,2,2)=mat3;
   %tk{k}(:,:,1,2)=mat3;
end
tt=tt_matrix(tt);
