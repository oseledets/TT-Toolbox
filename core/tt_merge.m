function [tt1]=tt_merge(tt,k)
%[TT1]=TT_MERGE(TT,K)
%Merges two dimensions K,K+1 in TT, result is a D-1 dimensional tensor
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
d=size(tt,1);
tt1=cell(d-1,1); %Dimension K increases, K+1 disappears
if ( k== 1 )
  mat1=tt{1};
  mat2=tt{2};
  n1=size(mat1,1);%r1=size(mat1,2);
  r3=size(mat2,3);
  n2=size(mat2,1);
  %mat2 is n2 x r1 x r3
  %mat1 is n1 x r1
  mat=ten_conv(mat2,2,mat1'); %mat is n2 x n1 x r3
  mat=permute(mat,[2,1,3]);
  mat=reshape(mat,n1*n2,r3);
  tt1{1}=mat;
tt1(k+1:d-1)=tt(k+2:d);
return
elseif( k > 1 )
 tt1(1:k-1)=tt(1:k-1);
else
  error('Incorrect K=%d, should be >= 1',k);    
end
core1=tt{k};
core2=tt{k+1};
n1=size(core1,1); r2=size(core1,2); r3=size(core1,3);
n2=size(core2,1); r4=size(core2,3);
core=reshape(core1,[n1*r2,r3])*reshape(permute(core2,[2,1,3]),[r3,n2*r4]);
%core is n1 x r2 x n2 x r4
core=permute(reshape(core,[n1,r2,n2,r4]),[1,3,2,4]);
tt1{k}=reshape(core,[n1*n2,r2,r4]);
tt1(k+1:d-1)=tt(k+2:d);

return
end
