function [res]=tt_hdm(tt1,tt2)
%[RES]=TT_HDM(TT1,TT2)
%Hadamard (pointwise) product of two TT-tensors
%Ranks are multiplied
%
%
% TT Toolbox 2.1, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%
%------------------------------
% Vectorized acceleration by dc
%---------------------------
  d=size(tt1,1);
  res=cell(d,1);
  mat1=tt1{1};
  mat2=tt2{1};
  ncur=size(mat1,1);
  r1=size(mat1,2);
  r3=size(mat2,2);
  mat1 = reshape(mat1, ncur*r1, 1);
  mat1 = mat1*ones(1,r3);
  mat1 = reshape(mat1, ncur, r1*r3);
  mat2 = reshape(mat2, ncur*r3, 1);
  mat2 = mat2*ones(1,r1);
  mat2 = permute(reshape(mat2, ncur, r3,r1), [1 3 2]);
  mat2 = reshape(mat2, ncur, r1*r3);  
  res{1} = mat1.*mat2;

  
  mat1=tt1{d};
  mat2=tt2{d};
  ncur=size(mat1,1);
  r1=size(mat1,2);
  r3=size(mat2,2);
  mat1 = reshape(mat1, ncur*r1, 1);
  mat1 = mat1*ones(1,r3);
  mat1 = reshape(mat1, ncur, r1*r3);
  mat2 = reshape(mat2, ncur*r3, 1);
  mat2 = mat2*ones(1,r1);
  mat2 = permute(reshape(mat2, ncur, r3,r1), [1 3 2]);
  mat2 = reshape(mat2, ncur, r1*r3);  
  res{d} = mat1.*mat2;
  
  for k=2:d-1
      core1=tt1{k}; core2=tt2{k}; ncur=size(core1,1); r1=size(core1,2);
      r2=size(core1,3); r3=size(core2,2); r4=size(core2,3);
      
      core1 = reshape(core1, ncur*r1*r2, 1);
      core1 = core1*ones(1,r3*r4);
      core1 = reshape(core1, ncur, r1*r2*r3*r4);
      
      core2 = reshape(core2, ncur*r3*r4, 1);
      core2 = core2*ones(1,r1*r2);
      core2 = reshape(core2, ncur, r3*r4, r1*r2);
      core2 = permute(core2, [1 3 2]);
      core2 = reshape(core2, ncur, r1*r2*r3*r4);
      
      res{k}=core1.*core2;
      res{k}=reshape(res{k}, ncur, r1, r2, r3, r4);
      res{k}=permute(res{k}, [1 2 4 3 5]);
      res{k}=reshape(res{k}, ncur, r1*r3, r2*r4);
  end
  return

end
