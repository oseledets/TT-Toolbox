function [res]=tt_hdm(tt1,tt2)
%[RES]=TT_HDM(TT1,TT2)
%Hadamard (pointwise) product of two TT-tensors
%Ranks are multiplied
%
%
% TT Toolbox 1.0, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
  d=size(tt1,1);
  res=cell(d,1);
  mat1=tt1{1};
  mat2=tt2{1};
  ncur=size(mat1,1);
  r1=size(mat1,2);
  r2=size(mat2,2);
  mat=zeros(ncur,r1,r2);
  for i=1:ncur
      for s1=1:r1
          for s2=1:r2
              mat(i,s1,s2)=mat1(i,s1)*mat2(i,s2);
          end
      end
  end
  res{1}=reshape(mat,[ncur,r1*r2]);
  mat1=tt1{d};
  mat2=tt2{d};
  ncur=size(mat1,1);
  r1=size(mat1,2);
  r2=size(mat2,2);
  mat=zeros(ncur,r1,r2);
  for i=1:ncur
      for s1=1:r1
          for s2=1:r2
              mat(i,s1,s2)=mat1(i,s1)*mat2(i,s2);
          end
      end
  end
  res{d}=reshape(mat,[ncur,r1*r2]);
  for k=2:d-1
      core1=tt1{k}; core2=tt2{k}; ncur=size(core1,1); r1=size(core1,2);
      r2=size(core1,3); r3=size(core2,2); r4=size(core2,3);
      core=zeros(ncur,r1,r3,r2,r4);
      for i=1:ncur
          for s1=1:r1
              for s2=1:r2
                  for s3=1:r3
                      for s4=1:r4
                          core(i,s1,s3,s2,s4)=core1(i,s1,s2)*core2(i,s3,s4);
                      end
                  end
              end
          end
      end
      res{k}=reshape(core,[ncur,r1*r3,r2*r4]);
  end
  return
