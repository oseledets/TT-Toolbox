function [y]=tt_crs2(arr,eps,nswp)
%[A]=TT_CRS2(ARR,EPS,NSWP)
%Cross approximation method for array ARR, 
%which tries to compute everything by
%increasing rank as far as possible
%In future will need the handler to the function
%NSWP - optional parameter, number of sweeps


%PARAMETERS SECTION


%If we increase rank at position i by 1 then it is easy to find 
%B(a1,i1,a2)*B(a2,i2,a3) -> ind1(a1,i1,a2), ind2(a2,i2,a3) 
%We compute Phi(a1,i1,i2,a3) and using ind1, ind2 we update them
%giving ind1(a1,i1,a2')*ind2(a2',ind2,a3) (two indices are updated)
%ind(a1,i1,a2) shows, how to select a2 rows out of (a1,i1) for the 
%approximation, so Schur complement approximation adds 3-mode indices
%to ind1 and second-mode ones to ind2, so this code is (completely) 
%incorrect (should be rewritten fully)

d=arr{1}; %Array size
sz=arr{2};

ind=cell(d,1); %This guy serves as (0,1) storage vector, showing how to 
%select from indices
for i=1:d
  ind{i}=zeros(1,sz(i),1);
  ind{i}(1) = 1;
end
rk=ones(1,d+1); %Ranks
%Forward sweep --- let us compute the solution in the "new" format
%initialized left_index & right_index 
%Compute right index sets for ind
ind0=ind{d}; ind0=reshape(ind0,[rk(d),sz(d)]);
ind_right=cell(d,1);
ind_right{d}=ind0*(1:sz(d))';
ind_left=cell(d,1);
ind0=ind{1}; ind0=reshape(ind0,[sz(1),rk(2)]);
ind_left{1}=(1:sz(1))*ind0;

for i=d-1:-1:2
   %Compute ind_right{i-1} from ind_right{i}
   %We need: r last part of ind_right{i}
   %and r first indices; how to get them from ind0 ??
   
   ind0=ind{i}; %ind0 is rk(i)*sz(i)*rk(i+1)
   %if: reshape into a (rk(i)*[sz(i)*rk(i+1)]) matrix, then 
   %for each r(k) we have a matrix n(i)*rk(i+1) with only one nonzero 
   %element. so, if we multiply it by say n(i)*ones.... we get required
   %n(i)
   %so we have: to convolve with 1:sz(i) over the second mode and
   %ones(1,rk(i+1)) over the third mode, whereas convolve over 
   rr=reshape(ind0,[rk(i)*sz(i),rk(i+1)])*((1:rk(i+1))'); 
   rr=reshape(rr,[rk(i),sz(i)]); rr=rr*((1:sz(i))');
   n=reshape(ind0,[rk(i)*sz(i),rk(i+1)])*(ones(rk(i+1),1)); 
   n=reshape(n,[rk(i),sz(i)]); n=n*ones(sz(i),1);
   ind_right{i}=[n,ind_right{i+1}(rr,:)]; %ind_right{i} has rk(i) elements
end

nswp = 5;
swp=1;
while ( swp < nswp )
   %Forward iteration
   %Compute matrix b
   %b(n1,n2,ra2)
   ind1=ind{1}; ind2=ind{2};
   n1=size(ind1,2); n2=size(ind2,2);
   ra2=size(ind2,3);
   ra1=1;
   b=zeros(n1,n2,ra2);
   for i1=1:n1
     for i2=1:n2
        for s=1:ra2
           ind_loc=[i1,i2,ind_right{3}];
           b(i1,i2,s)=my_elem(ind_loc,arr);
        end
     end
   end
   %ind1 and ind2 constitutes cross approximation to b; 
   b=reshape(b,[ra1*n1,n2*ra2]); 
   ind1=reshape(ind1,[n1,rk(2)]);
   ind1=(1:n1)*ind1;
   ind2=reshape(ind2,[rk(2),n2*ra2]);
   ind2=ind2*((1:n2*ra2)');
   %Now we can extend ind1 & ind2 by Schur complement technique   
   ind_lf1=setdiff(1:ra1*n1,ind1);
   ind_rf1=setdiff(1:n2*ra2,ind2);

          sbm=b(ind1,ind2);
     
       b_schr=b(ind_lf1,ind_rf1)-(b(ind_lf1,ind2)/sbm)*b(ind1,ind_rf1);
       %Now approximate me by skeleton method --- only surrogate now
       nrm=norm(b,'fro');
       [u,s,v]=svd(b_schr,'econ');
       s=diag(s); 
       
       if ( norm(b_schr,'fro') > eps*nrm )
       
          r=my_chop2(s,eps*nrm);
          u=u(:,1:r); v=v(:,1:r);
          ind_lf2=maxvol2(u); ind_rf2=maxvol2(v);
          ind1=[ind1,ind_lf1(ind_lf2)];
          ind2=[ind2,ind_rf1(ind_rf2)];
       end
       %Convert ind1 & ind2 to "old" ind format
       rnew=size(ind1,2); 
       ind1_new=zeros(ra1*n1,rnew); 
       ind2_new=zeros(rnew*n2,ra2);
       ind2_old=ind2{2};
       for s=1:rnew
           ind1_new(ind1(s),s)=1;
       end
       ind{1}=reshape(ind1_new,[ra1,n1,rnew]);
       ind{2}=reshape(ind2_new,[rnew,n2,ra2]);
       ind_left{1}=update_index_left(1,ind{1},[]);
       ind_left{2}=update_index_left(2,ind{2},ind_left{1});
       keyboard;
       rk(2)=rnew;
       %Now compute NEW ind_left ind_left{i-1} & ind_left{i} are modified
       %:( 
       i=2;
       ind0=ind{i-1}; %ind0 is rk(i)*sz(i)*rk(i+1)
      %if: reshape into a (rk(i)*[sz(i)*rk(i+1)]) matrix, then 
      %for each r(k) we have a matrix n(i)*rk(i+1) with only one nonzero 
      %element. so, if we multiply it by say n(i)*ones.... we get required
      %n(i)
   %so we have: to convolve with 1:sz(i) over the second mode and
   %ones(1,rk(i+1)) over the third mode, whereas convolve over 
   rr=(1:rk(i-1))*reshape(ind0,[rk(i-1),sz(i)*rk(i)]); 
   rr=reshape(rr,[sz(i),rk(i)]); rr=(1:sz(i))*rr;
   n=(1:rk(i-1))*reshape(ind0,[rk(i-1),sz(i)*rk(i)]); 
   n=reshape(n,[sz(i),rk(i)]); n=ones(1,sz(i))*n;
   ind_left{i}=[ind_left{i-1}(rr,:),n]; %ind_left{i} has rk(i+1) elements 
      keyboard;
   for i=2:d-1
       ind1=ind{1}; ind2=ind{2};
   n1=size(ind1,2); n2=size(ind2,2);
   ra2=size(ind2,3);
   ra1=1;
   b=zeros(n1,n2,ra2);
   for i1=1:n1
     for i2=1:n2
        for s=1:ra2
           ind_loc=[i1,i2,ind_right{3}];
           b(i1,i2,s)=my_elem(ind_loc,arr);
        end
     end
   end
   %ind1 and ind2 constitutes cross approximation to b; 
   b=reshape(b,[ra1*n1,n2*ra2]); 
   ind1=reshape(ind1,[n1,rk(2)]);
   ind1=(1:n1)*ind1;
   ind2=reshape(ind2,[rk(2),n2*ra2]);
   ind2=ind2*((1:n2*ra2)');
   %Now we can extend ind1 & ind2 by Schur complement technique   
   ind_lf1=setdiff(1:ra1*n1,ind1);
   ind_rf1=setdiff(1:n2*ra2,ind2);

          sbm=b(ind1,ind2);
     
       b_schr=b(ind_lf1,ind_rf1)-(b(ind_lf1,ind2)/sbm)*b(ind1,ind_rf1);
       %Now approximate me by skeleton method --- only surrogate now
       nrm=norm(b,'fro');
       [u,s,v]=svd(b_schr,'econ');
       s=diag(s); 
       
       if ( norm(b_schr,'fro') > eps*nrm )
       
          r=my_chop2(s,eps*nrm);
          u=u(:,1:r); v=v(:,1:r);
          ind_lf2=maxvol2(u); ind_rf2=maxvol2(v);
          ind1=[ind1,ind_lf1(ind_lf2)];
          ind2=[ind2,ind_rf1(ind_rf2)];
       end
       %Convert ind1 & ind2 to "old" ind format
       rnew=size(ind1,2); 
       ind1_new=zeros(ra1*n1,rnew); 
       ind2_new=zeros(rnew,n2*ra2);
       for s=1:rnew
           ind1_new(ind1(s),s)=1;
           ind2_new(s,ind2(s))=1;
       end
       ind{1}=reshape(ind1_new,[ra1,n1,rnew]);
       ind{2}=reshape(ind2_new,[rnew,n2,ra2]);
       
   end
       %Now construct full approximation
       %sbm=b(ind1,ind2);
       %u=(b(:,ind2)/sbm);
       %v=b(ind1,:);
   swp=swp+1;
end
return
end

function [new_ind]=update_index_left(num,ind,old_ind)
  if ( num == 1 )
      n1=size(ind,2); rk=size(ind,3);
     ind=reshape(ind,[n1,rk]);
    new_ind=(1:n1)*ind;

  else
      n1=size(ind,2); rk1=size(ind,1); rk2=size(ind,3);
     rr=(1:rk1)*reshape(ind,[rk1*n1,rk2]); 
   rr=reshape(rr,[n1,rk2]); rr=(1:n1)*rr;
   n=(1:rk1)*reshape(ind,[rk1,n1*rk2]); 
   n=reshape(n,[n1,rk2]); n=ones(1,n1)*n;
   new_ind=[old_ind(rr,:),n]; %ind_left{i} has rk(i+1) elements 
 
  end
return
end