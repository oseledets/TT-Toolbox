function [res]=tt_to_qtt(tt,eps,d0)
%[RES]=TT_TO_QTT(TT,EPS,D0)
% Converts TT format to QTT
% RES is the QTT representation for the core
% EPS is the allowed accuracy
%
%
%
% TT Toolbox 1.0, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

%d0 is the QTT dimension;
%Algorithm.
%1) Compute left-to-right QR
%2) Start computing right-to-left QR, 
%additionally compressing obtained cores to QTT
d=size(tt,1);
res=cell(d*d0,1); %Too much?
tt1=tt;
%Compute right-to-left qr
mat=tt1{d};
[q,r]=qr(mat,0);
%r_right{d}=r;
tt1{d}=q;
for i=(d-1):-1:2
    core=tt1{i};
    core=ten_conv(core,3,r');
    ncur=size(core,1);
    r2=size(core,2);
    r3=size(core,3);
    core=permute(core,[1,3,2]);
    core=reshape(core,[ncur*r3,r2]);
    [tt1{i},r]=qr(core,0);
    %r_right{i}=r;
    rnew=min(r2,ncur*r3);
    tt1{i}=reshape(tt1{i},[ncur,r3,rnew]);
    tt1{i}=permute(tt1{i},[1,3,2]);
end
tt1{1}=tt1{1}*r';
mat=tt1{1};
%Compress mat to QTT
r2=size(mat,2);
mat=reshape(mat,[2*ones(1,d0),r2]);
%tts=full_to_tt(mat,eps/sqrt(d-1));
tts=full_to_tt(mat,-1.0);
res(1:d0)=tts(1:d0);
r=tts{d0+1}';
%[q,r]=qr(mat,0);
%r_right{d}=r;
%tt1{d}=q;
for i=2:d-1
    core=tt1{i};
    core=ten_conv(core,2,r');
    ncur=size(core,1);
    r2=size(core,2);
    r3=size(core,3);
    core=permute(core,[2,1,3]);
    core=reshape(core,[r2,2*ones(1,d0),r3]);
    %tts=full_to_tt(core,eps/sqrt(d-1));
    tts=full_to_tt(core,-1.0);
    res(((i-1)*d0+1):i*d0)=tts(2:(d0+1));
    %keyboard;
    res{(i-1)*d0+1}=ten_conv(res{(i-1)*d0+1},2,tts{1}');
    r=tts{d0+2}';
    %core=reshape(core,[ncur*r2,r2]);
    %[tt1{i},r]=qr(core,0);
    %r_right{i}=r;
    %rnew=min(r2,ncur*r3);
    %tt1{i}=reshape(tt1{i},[ncur,r3,rnew]);
    %tt1{i}=permute(tt1{i},[1,3,2]);
end
mat=tt1{d}*r';
r2=size(mat,2); 
mat=mat';
mat=reshape(mat,[r2,2*ones(1,d0)]);
tts=full_to_tt(mat,-1.0);%eps1/sqrt(d-1));
res(((d-1)*d0+1):(d*d0))=tts(2:(d0+1));
res{(d-1)*d0}=ten_conv(res{(d-1)*d0},3,tts{1}');
%keyboard;
%[factors{1},s,res{1}]=svd(tt{1},'econ');
%res{1}=s*res{1}';
%[factors{d},s,res{d}]=svd(tt{d},'econ');
%res{d}=s*res{d}';

%Compute reduced cores & their QTT decomposition
for i=2:d-1
  core=tt{i};
  core=ten_conv(core,2,r_left{i-1}');
  core=ten_conv(core,3,r_right{i+1}');
  ncur=size(core,1); r2=size(core,2); r3=size(core,3);
    %keyboard;
  core1=core;
  core1=reshape(core1,[2*ones(1,10),r2,r3]);
  %tt0=full_to_tt(core1,eps);
  %fprintf('i=%d, rank=%3.1f \n',i,tt_erank(tt0));
  %keyboard;
  core=reshape(core,[ncur,r2*r3]);
  [u,s,v]=svd(core,'econ');
  r=my_chop(diag(s),eps/sqrt(d-1));
 u=u(:,1:r); 
 factors{i}=u;
  %Compute cores as convolution
  res{i}=ten_conv(tt{i},1,u);
  %res{i}=s(1:r,1:r)*v(:,1:r)'; 
  %res{i}=reshape(res{i},r,r2,r3);
end
return
end
