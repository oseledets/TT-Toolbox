function [tt,ind_right]=tt_crossr(arr,ind_left)
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

d=arr{1};
sz=arr{2};
%arr1=arr(:);
tt=cell(d,1);
ind_right=cell(d,1); %Computed ind_right
ind_l=ind_left{d-1};
r1=size(ind_l,2);
ncur=sz(1);
mat=zeros(ncur,r1);
for j=1:ncur
   for s=1:r1
      ind_f=[ind_l(:,s)',j];
      %mat(j,s)=arr1(sub_to_ind(ind_f,sz));
      val=my_elem(ind_f,arr); 
      mat(j,s)=val;
   end
end
[qs,rs]=qr(mat,0);
ind=maxvol2(qs);
ind_right{d-1}=ind;
%mat=mat*inv(mat(ind,:));
mat=qs*inv(qs(ind,:));
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
                 val=my_elem(ind_f,arr); 
                 core(i,s2,s3)=val;
                 %core(i,s2,s3)=arr1(sub_to_ind(ind_f,sz));

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
     %core=core*inv(core(ind,:)); 
     core=qs*inv(qs(ind,:));
     
     % core=reshape(core,[r2,ncur,r3]); core=permute(core,[2,1,3]);
     % core=reshape(core,[ncur,r3,r2]); core=permute(core,[1,3,2]);
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
      val=my_elem(ind_f,arr);
      mat(i,s)=val;
      %mat(i,s)=arr1(sub_to_ind(ind_f,sz));
    end
end
tt{1}=mat;
%for k=1:d
%     res=tt{k}-tt_true{k};
%    res=res(:);
%    fprintf('k=%d, approximation error: %e \n',k,norm(res));
%end
return
end
