function [res]=tt_kron2(tt1,tt2)
%[RES]=TT_KRON2(TT1,TT2)
% Computes Kronecker product between TT1 & TT2
% Cores are stored in a "new" enumeration
% Kronecker products unites indices as 
% (i1,j1) (i2,j2) and the mode sizes are multiplied
% Number of dimensions in tt1 and tt2 should be equal!
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
d1=size(tt1,1);
d2=size(tt2,1);
if ( d1 ~= d2 )
  fprintf('Size of tt1 is not equal to the size tt2 \n');
  return
end
d=d1;
res=cell(d,1);
mat1=tt1{1}; mat2=tt2{1}; n=size(mat1,1); rn=size(mat1,2);
m=size(mat2,1); rm=size(mat2,2);
mat=kron(reshape(mat1,[n*rn,1]),reshape(mat2,[m*rm,1]));
mat=reshape(mat,[n,rn,m,rm]); mat=permute(mat,[1,3,2,4]);
mat=reshape(mat,[n*m,rn*rm]);
res{1}=mat;
mat1=tt1{d}; mat2=tt2{d}; n=size(mat1,2); rn=size(mat1,1);
m=size(mat2,2); rm=size(mat2,1);
mat=kron(reshape(mat1,[1,rn*n]),reshape(mat2,[1,rm*m]));
mat=reshape(mat,[rn,n,rm,m]); mat=permute(mat,[1,3,2,4]);
mat=reshape(mat,[rn*rm,n*m]);
res{d}=mat;
for k=2:d-1
  core1=tt1{k}; core2=tt2{k};
  n=size(core1,2); m=size(core2,2);
  rn1=size(core1,1); rn2=size(core1,3);
  rm1=size(core2,1); rm2=size(core2,3);
  core=kron(reshape(core1,[rn1*n*rn2,1]),reshape(core2,[rm1*m*rm2,1]));
  core=reshape(core,[rn1,n,rn2,rm1,m,rm2]);
  core=permute(core,[1,4,3,5,2,6]);
  core=reshape(core,[rn1*rm1,n*m,rn2*rm2]);
  res{k}=core;
end
return
end
