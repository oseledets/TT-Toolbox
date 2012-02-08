function [factors,res]=tt_tuck(tt,eps)
%Computes Tucker factors and Tucker core of the TT-tensor
%   [FACTORS,RES]=TT_TUCK(TT) Compute Tucker factors of the TT tensor 
%   FACTORS is the cell array of Tucker factors RES is the TT 
%   representation for the core, EPS is the allowed accuracy
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
tt=core(tt);
d=size(tt,1);

tt1=tt;
mat=tt1{1};
[q,r]=qr(mat,0);
tt1{1}=q; 
r_left=cell(d,1);
r_left{1}=r;
for i=2:d-1
    core1=tt1{i};
    core1=ten_conv(core1,2,r');
    ncur=size(core1,1);
    r2=size(core1,2);
    r3=size(core1,3);
    core1=reshape(core1,[ncur*r2,r3]);
    [tt1{i},r]=qr(core1,0);
    r_left{i}=r;
    rnew=min(ncur*r2,r3);
    tt1{i}=reshape(tt1{i},[ncur,r2,rnew]);
end
%Compute right-to-left qr & maxvol
tt1=tt;
mat=tt1{d};
r_right=cell(d,1);
[q,r]=qr(mat,0);
r_right{d}=r;
tt1{d}=q;
for i=(d-1):-1:2
    core1=tt1{i};
    core1=ten_conv(core1,3,r');
    ncur=size(core1,1);
    r2=size(core1,2);
    r3=size(core1,3);
    core1=permute(core1,[1,3,2]);
    core1=reshape(core1,[ncur*r3,r2]);
    [tt1{i},r]=qr(core1,0);
    r_right{i}=r;
    rnew=min(r2,ncur*r3);
    tt1{i}=reshape(tt1{i},[ncur,r3,rnew]);
    tt1{i}=permute(tt1{i},[1,3,2]);
end
%keyboard;
factors=cell(d,1);
res=cell(d,1);
[factors{1},s,res{1}]=svd(tt{1},'econ');
res{1}=s*res{1}';
[factors{d},s,res{d}]=svd(tt{d},'econ');
res{d}=s*res{d}';
%Compute reduced cores & Tucker factors
for i=2:d-1
  core1=tt{i};
  core1=ten_conv(core1,2,r_left{i-1}');
  core1=ten_conv(core1,3,r_right{i+1}');
  ncur=size(core1,1); r2=size(core1,2); r3=size(core1,3);
  core1=reshape(core1,[ncur,r2*r3]);
  [u,s,v]=svd(core1,'econ');
  r=my_chop2(diag(s),eps/sqrt(d-1)*norm(diag(s)));
 u=u(:,1:r); 
 factors{i}=u;
  %Compute cores as convolution
  res{i}=ten_conv(tt{i},1,u);
  %res{i}=s(1:r,1:r)*v(:,1:r)'; 
  %res{i}=reshape(res{i},r,r2,r3);
end
res=tt_tensor(res);
return
end
