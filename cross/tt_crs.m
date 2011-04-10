function [y]=tt_crs(arr,eps,nswp)
%[A]=TT_CRS(ARR,EPS,NSWP)
%Cross approximation method for array ARR, 
%which tries to compute everything by
%increasing rank as far as possible
%In future will need the handler to the function
%NSWP - optional parameter, number of sweeps
%PARAMETERS SECTION
rmin=1; verb=false;
radd=0;
if ( nargin == 2 )
nswp=2;
end
yold=[];
d=arr{1}; %Array size
sz=arr{2};

%verb=true;
y=cell(d,1);
%Init solution to start from some nonzero element of A
%Stupid initialization is to set ind = [1,1,1,1,]
start_ind=ones(1,d);
%And compute corresponding crosses
for i=1:d
   y{i}=zeros(sz(i),1);
   for k=1:sz(i)
      ind = [start_ind(1:i-1),k,start_ind(i+1:d)];
      y{i}(k)=my_elem(ind,arr);
   end
end
%y=tt_random(sz,d,4);
%Also here we need to evaluate "right-to-left" submatrices of y and
%initialized y_right

%Here we need: "left indices" and right indices; We consequently 
%orthogonalize y from left to right (assuming that from right-to-left 
%they are already orthogonalized)
%At this point we assume that: Y is orthogonal from right-to-left
%and ind_right correspond to correct max-vol positions in this matrix
ind_left=cell(d,1);
er=2*eps;
swp=1;
while ( swp < nswp && er > eps )
%Compute right-to-left qr & maxvol
mat=y{d};
[q,rv]=qr(mat,0);
y{d}=q;
for i=(d-1):-1:2
    core=y{i};
    core=ten_conv(core,3,rv');
    ncur=size(core,1);
    r2=size(core,2);
    r3=size(core,3);
    core=permute(core,[1,3,2]);
    core=reshape(core,[ncur*r3,r2]);
    [y{i},rv]=qr(core,0);
    rnew=min(r2,ncur*r3);
    y{i}=reshape(y{i},[ncur,r3,rnew]);
    y{i}=permute(y{i},[1,3,2]);
end
y{1}=y{1}*rv';
sbm=cell(d,1);
mat=y{d};
ind_right=cell(d,1);
[ind]=maxvol2(mat);
ind_right{d-1}=ind;
r1=mat(ind,:);
ind_r=cell(d,1);
y{d}=mat/r1;
ind_r{d-1}=ind;
sbm{d}=r1;
for i=(d-1):-1:2
    core=y{i};
    ncur=size(core,1);
    r2=size(core,2);
    r3=size(core,3);
    core=ten_conv(core,3,r1');
    core=permute(core,[3,1,2]);
    core=reshape(core,[r3*ncur,r2]);
    [ind]=maxvol2(core); 
    ind_old=ind_right{i};
    rnew=min(ncur*r3,r2);
    ind_new=zeros(d-i+1,rnew);
     for s=1:rnew
        f_in=ind(s);
        w1=tt_ind2sub([r3,ncur],f_in);
        rs=w1(1); js=w1(2);
        ind_new(:,s)=[js,ind_old(:,rs)'];
     end
     ind_right{i-1}=ind_new;
    r1=core(ind,:);
    ind_r{i-1}=ind;
    sbm{i}=r1;
    y{i}=core/r1;
    rnew=min(ncur*r3,r2);
    y{i}=reshape(y{i},[r3,ncur,rnew]);
    y{i}=permute(y{i},[2,3,1]);
end

    
    
    
    %First core: we compute B(i1,i2,ra)
    n1=sz(1);
    n2=sz(2);
    ra2=size(ind_right{2},2);
    ind2=ind_right{2};
    b=zeros(n1,n2,ra2);
    for i=1:n1
      for j=1:n2
         for s2=1:ra2
             ind=[i,j,ind2(:,s2)'];
            b(i,j,s2)=my_elem(ind,arr);
         end
      end
    end
    %b=reshape(b,[n1*n2,ra]);
    %b=b/inv(sbm{3});
    %Compute its SVD, new rank estimate, and new first core
    b=reshape(b,[n1,n2*ra2]);
    
    [u,s,v]=svd(b,'econ');
    s=diag(s); nrm=norm(s);
    rold=size(y{1},2);
    r=my_chop2(s,eps*nrm);
    r=max([rmin,r,rold+1]);
    r=min(r,size(s,1));
    if ( verb ) 
        fprintf('We can push rank %d to %d \n',1,r);  
    end
    u=u(:,1:r);
    ind=maxvol2(u);
    ind_left{1}=ind;
    u=u/u(ind,:);
    y{1}=u;

    %Now check all other cores
   for i=2:d-2
       %We compute B(ra1,n1,n2,ra2)
       n1=sz(i); n2=sz(i+1);
       ra1=size(ind_left{i-1},2);
       ra2=size(ind_right{i+1},2);
       ind1=ind_left{i-1};
       ind2=ind_right{i+1};
       b=zeros(ra1,n1,n2,ra2);
       for i1=1:n1
         for i2=1:n2
           for s1=1:ra1
             for s2=1:ra2
                ind=[ind1(:,s1)',i1,i2,ind2(:,s2)'];
                b(s1,i1,i2,s2)=my_elem(ind,arr);
             end
           end
         end
       end
       b=reshape(b,[ra1*n1,n2*ra2]);
       [u,s,v]=svd(b,'econ'); s=diag(s);
      rold=size(y{i},3);
      r=my_chop2(s,eps*norm(s));
      r=max([r,rmin,rold+1]);
        r=min(r,size(s,1));
  %if ( 
  if ( verb )
    fprintf('We can push rank %d to %d \n',i,r);
  end
  u=u(:,1:r);
  %v=v(:,1:r)*diag(s(1:r));
  u=reshape(u,[ra1*n1,r]);
  rnew=min(ra1*n1,r);
  [u,rm]=qr(u,0);
  ind=maxvol2(u);
  ind_old=ind_left{i-1};
  ind_new=zeros(i,rnew);
    for s=1:rnew
     f_in=ind(s);
     w1=tt_ind2sub([ra1,n1],f_in);
     rs=w1(1); js=w1(2);
     ind_new(:,s)=[ind_old(:,rs)',js];
    end
    ind_left{i}=ind_new;
    %core=core*inv(core(ind,:)); 
    u=u/u(ind,:);
    u=reshape(u,[ra1,n1,rnew]); u=permute(u,[2,1,3]);
    y{i}=u;
   end
    %Now compute the last core - and we are done :-)
    %Compute B(ra1,n1,n2)
    n1=sz(d-1);
    n2=sz(d);
    ra1=size(ind_left{d-2},2);
    ind1=ind_left{d-2};
    b=zeros(ra1,n1,n2);
    for i1=1:n1
      for i2=1:n2
         for s1=1:ra1
            ind=[ind1(:,s1)',i1,i2];
            b(s1,i1,i2)=my_elem(ind,arr);
         end
      end
    end
    b=reshape(b,[ra1*n1,n2]);
    [u,s,v]=svd(b,'econ'); s=diag(s);
    
    rold=size(y{i},3);
    r=my_chop2(s,eps*norm(s));
    r=max([r,rmin,rold+radd]);
    r=min(r,size(s,1));
    if ( verb )
      fprintf('We can push rank %d to %d \n',d-1,r);
    end
    u=u(:,1:r);
    v=v(:,1:r);
    s=s(1:r);
    u=reshape(u,[ra1,n1,r]);
    u=permute(u,[2,1,3]);
    y{d-1}=u;
    y{d}=v*diag(s);
    if ( ~isempty(yold) ) 
      er = tt_dist2(yold,y)/sqrt(tt_dot(y,y));
      fprintf('er=%3.2e \n',er);
    end
    yold=y;
   swp=swp+1;
end

return
end
