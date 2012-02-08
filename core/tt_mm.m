function [res]=tt_mm(mat1,mat2)
%Matrix-by-matrix product in TT1.0 format
%   [RES]=TT_MM(MAT1,MAT2) Matrix MAT1 by matrix MAT2 product in the 
%   TT-format. Please avoid its usage: it will be removed in
%   future releases. Use * operator from the object-oriented version.
%
%
% TT-Toolbox 2.2, 2009-2011
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
d=size(mat1,1);
res=cell(d,1);
k=1;
for s=1:2
n1=size(mat1{k},1);
m1=size(mat1{k},2);
rm=size(mat1{k},3);
n2=size(mat2{k},1);
m2=size(mat2{k},2);
rv=size(mat2{k},3);
%mat1{1} is (n1 x m1 x rm , mat2{1} is n2xm2xrm, result is n1xm2xrv1*rv)
%res{1}=reshape(permute(mat1{1},[1,3,2],[n1*rm1]),m1)*reshape(mat2{1},[n2,m2*rv2]);
res{k}=reshape(permute(mat1{k},[1,3,2]),[n1*rm,m1])*reshape(mat2{k},[n2,m2*rv]);
%res{1} is n1 x rm x m2 x rv
res{k}=reshape(res{k},[n1,rm,m2,rv]);
res{k}=permute(res{k},[1,3,2,4]);
res{k}=reshape(res{k},[n1,m2,rm*rv]);
if ( k == 1 ) 
 k=d;
end
end
for i=2:d-1
  n1=size(mat1{i},1);
  m1=size(mat1{i},2);
  rm1=size(mat1{i},3);
  rm2=size(mat1{i},4);
  n2=size(mat2{i},1);
  m2=size(mat2{i},2);
  rv1=size(mat2{i},3);
  rv2=size(mat2{i},4);
  %mat1{i} is n1 x m1 x rm1 x rm2 
  %mat2{i} is n2 x m2 x rv1 x rv2
  %n1 x rm1 x rm2 x m1 * n2 x m2 x rv1 x rv2
  %n1 x rm1 x rm2 x m2 x rv1 x rv2
  %Need: n1 x m2 x rm1 x rv1 x rm2 x rv2 
  %Permut: 1 4 2 5 3 6
  res{i}=reshape(permute(mat1{i},[1,3,4,2]),[n1*rm1*rm2,m1])*reshape(mat2{i},[n2,m2*rv1*rv2]);
  res{i}=reshape(res{i},[n1,rm1,rm2,m2,rv1,rv2]);
  res{i}=permute(res{i},[1,4,2,5,3,6]);
  res{i}=reshape(res{i},[n1,m2,rm1*rv1,rm2*rv2]);
end 
return
end
