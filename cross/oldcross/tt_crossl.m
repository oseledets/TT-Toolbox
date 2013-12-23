function [tt,ind_left]=tt_crossl(d,n,fun,ind_right)
%One left-to-right sweep of the TT-cross method.
%   [TT,IND_LEFT]=TT_CROSSL(D,N,FUN,IND_RIGHT) Computes one QR-maxvol sweep
%   of the TT-cross method. The input is the pair (D,N) that determines the
%   size of the tensor, FUN is the function handle to evaluate a particular
%   element of the tensor, i.e., VAL=FUN(IND). To pass parameters, please
%   use anonymous function handles in MATLAB
% 
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

if ( numel(n) == 1 )
  n=n*ones(d,1);
end
sz=n;
tt=cell(d,1);
ind_left=cell(d,1); %Computed ind_left
ind_r=ind_right{1};
r1=size(ind_r,2);
ncur=sz(1);
mat=zeros(ncur,r1);
for i=1:ncur
   for s=1:r1
      ind_f=[i,ind_r(:,s)'];
      val=fun(ind_f);
      mat(i,s)=val;
   end
end
[qs,~]=qr(mat,0);
ind=maxvol2(qs);
ind_left{1}=ind;
mat=qs/qs(ind,:);
tt{1}=mat;

for k=2:d-1
     ncur=sz(k);
     ind_l=ind_left{k-1};
     ind_r=ind_right{k};
     r2=size(ind_l,2);
     r3=size(ind_r,2);
     core=zeros(ncur,r2,r3);
     for i=1:ncur
         for s2=1:r2
             for s3=1:r3
                 ind_f=[ind_l(:,s2)',i,ind_r(:,s3)'];
                 val=fun(ind_f);
                 core(i,s2,s3)=val;
             end
         end
     end
     core=permute(core,[2,1,3]); core=reshape(core,[r2*ncur,r3]);
     rnew=min(r2*ncur,r3);
     [qs,rs]=qr(core,0);
     ind=maxvol2(qs);
     ind_old=ind_left{k-1};
     ind_new=zeros(k,rnew);
     ncur=sz(k);
     for s=1:rnew
        f_in=ind(s);
        w1=tt_ind2sub([r2,ncur],f_in);
        rs=w1(1); js=w1(2);
        ind_new(:,s)=[ind_old(:,rs)',js];
     end
     ind_left{k}=ind_new;
     core=qs/qs(ind,:);
     core=reshape(core,[r2,ncur,rnew]); core=permute(core,[2,1,3]);
   
    tt{k}=core;
end
ncur=sz(d);
ind_l=ind_left{d-1};
r=size(ind_l,2);
mat=zeros(ncur,r);
for j=1:ncur
    for s=1:r
      ind_f=[ind_l(:,s)',j];
      val=fun(ind_f);
      mat(j,s)=val;
    end
end
tt{d}=mat;
tt=tt_tensor(tt);
return
end
