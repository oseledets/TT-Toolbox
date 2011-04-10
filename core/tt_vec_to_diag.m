function [tt_mat]=tt_vec_to_diag(tt_vec)
%Converts tt-vector to a diagonal matrix in tt-format.
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
d=size(tt_vec,1);
tt_mat=cell(d,1);
r=size(tt_vec{1},2);
n=size(tt_vec{1},1);
tt_mat{1}=zeros(n,n,r);
for i=1:r
    tt_mat{1}(:,:,i)=diag(tt_vec{1}(:,i));
end

r=size(tt_vec{d},2);
n=size(tt_vec{d},1);
tt_mat{d}=zeros(n,n,r);
for i=1:r
    tt_mat{d}(:,:,i)=diag(tt_vec{d}(:,i));
end

for i=2:d-1
   r2=size(tt_vec{i},2);
   r3=size(tt_vec{i},3);
    n=size(tt_vec{i},1);
    tt_mat{i}=zeros(n,n,r2,r3);
    for j=1:r2
        for l=1:r3
            tt_mat{i}(:,:,j,l)=diag(tt_vec{i}(:,j,l));
        end
    end
end

return
end
