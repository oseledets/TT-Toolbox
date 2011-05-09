%Seconds technical function that tries to mimic the skeleton decomposition
%of the TT-tensor directly using properly computed index sets
function [tt,ind_left]=tt_crossl(arr,ind_right)
d=arr{1};
sz=arr{2};
tt=cell(d,1);
ind_left=cell(d,1); %Computed ind_left
ind_r=ind_right{1};
r1=size(ind_r,2);
ncur=sz(1);
mat=zeros(ncur,r1);
for i=1:ncur
   for s=1:r1
      ind_f=[i,ind_r(:,s)'];
      val=my_elem(ind_f,arr);
      mat(i,s)=val;
   end
end
[qs,rs]=qr(mat,0);
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
                 val=my_elem(ind_f,arr);
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
     %core=core*inv(core(ind,:)); 
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
      val=my_elem(ind_f,arr);
      mat(j,s)=val;
      %mat(j,s)=arr1(sub_to_ind(ind_f,sz));
    end
end
tt{d}=mat;
%for k=1:d
%     res=tt{k}-tt_true{k};
%    res=res(:);
%    fprintf('k=%d, approximation error: %e \n',k,norm(res));
%end
return
end
