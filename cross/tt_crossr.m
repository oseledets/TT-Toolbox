function [tt,ind_right]=tt_crossr(d,n,fun,ind_left)
%One right-to-left sweep of the TT-cross method.
%   [TT,IND_RIGHT]=TT_CROSSR(D,N,FUN,IND_LEFT) Computes one QR-maxvol sweep
%   of the TT-cross method. The input is the pair (D,N) that determines the
%   size of the tensor, FUN is the function handle to evaluate a particular
%   element of the tensor, i.e., VAL=FUN(IND). To pass parameters, please
%   use anonymous function handles in MATLAB. 
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
ind_right=cell(d,1); %Computed ind_right
ind_l=ind_left{d-1};
r1=size(ind_l,2);
ncur=sz(1);
mat=zeros(ncur,r1);
for j=1:ncur
   for s=1:r1
      ind_f=[ind_l(:,s)',j];
      val=fun(ind_f); 
      mat(j,s)=val;
   end
end
[qs,rs]=qr(mat,0);
ind=maxvol2(qs);
ind_right{d-1}=ind;
%mat=mat*inv(mat(ind,:));
mat=qs/qs(ind,:);
tt{d}=mat;

for k=(d-1):-1:2
     ncur=sz(k);
     ind_l=ind_left{k-1};
     ind_r=ind_right{k};
     r2=size(ind_l,2);
     r3=size(ind_r,2);
     core=zeros(ncur,r2,r3);
     %fprintf('k=%d \n');
     for i=1:ncur
         for s2=1:r2
             for s3=1:r3
                 ind_f=[ind_l(:,s2)',i,ind_r(:,s3)'];
                 val=fun(ind_f); 
                 core(i,s2,s3)=val;

             end
         end
     end
     %core is ncurxr2xr3, from the left we need 
     %r2xncurxr3, but from the right 
     %we need (ncurxr3)xr2 or (r3xncur)xr2
     core=permute(core,[3,1,2]); core=reshape(core,[r3*ncur,r2]);
     rnew=min(r2,r3*ncur);
     [qs,rs]=qr(core,0);
     ind=maxvol2(qs);
     ind_old=ind_right{k};
     ind_new=zeros(d-k+1,rnew);
     ncur=sz(k);
     for s=1:rnew
        f_in=ind(s);
        w1=tt_ind2sub([r3,ncur],f_in);
        rs=w1(1); js=w1(2);
        ind_new(:,s)=[js,ind_old(:,rs)'];
     end
     ind_right{k-1}=ind_new;
     core=qs/(qs(ind,:));
     
     core=reshape(core,[r3,ncur,rnew]); core=permute(core,[2,3,1]);
   
    tt{k}=core;
end
ncur=sz(1);
ind_r=ind_right{1};
r=size(ind_r,2);
mat=zeros(ncur,r);
for i=1:ncur
    for s=1:r
      ind_f=[i,ind_r(:,s)'];
      val=fun(ind_f);
      mat(i,s)=val;
    end
end
tt{1}=mat;

tt=tt_tensor(tt);
return
end
