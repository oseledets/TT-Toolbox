function [y]=tt_crs4(arr,eps,nswp,y)
%[A]=TT_CRS4(ARR,EPS,NSWP)
%Cross approximation method for array ARR, 
%which tries to compute everything by
%increasing rank as far as possible
%In future will need the handler to the function
%NSWP - optional parameter, number of sweeps

%PARAMETERS SECTION
rmin=1; 
verb=true;
radd=0;
if ( nargin == 2 )
nswp=10;
end
yold=[];
d=arr{1}; %Array size
sz=arr{2};
%Compute full array
  %full_arr=zeros(sz);
  %for ind=1:prod(sz)
  %   id=tt_ind2sub(sz,ind);
  %   full_arr(ind) = my_elem(id,arr);
  %end
  %fa=full_arr(:);
  %B=arr{6};
  %keyboard;
%verb=true;
%y=cell(d,1);
%Init solution to start from some nonzero element of A
%Stupid initialization is to set ind = [1,1,1,1,]
%start_ind=1*ones(1,d);
%And compute corresponding crosses
%for i=1:d
%   y{i}=zeros(sz(i),1);
%   for k=1:sz(i)
%      ind = [start_ind(1:i-1),k,start_ind(i+1:d)];
%      y{i}(k)=my_elem(ind,arr);
%   end
%end
if (nargin < 4 )
y=tt_random(sz,d,5); y=tt_compr2(y,1e-8);
end

%y=tt_ones(d,sz);
%y=tt_compr2(y,1e-8);
%Also here we need to evaluate "right-to-left" submatrices of y and
%initialized y_right

%Here we need: "left indices" and right indices; We consequently 
%orthogonalize y from left to right (assuming that from right-to-left 
%they are already orthogonalized)
%At this point we assume that: Y is orthogonal from right-to-left
%and ind_right correspond to correct max-vol positions in this matrix

%Warmup procedure
ind_left=cell(d,1);
er=2*eps;
swp=1;

sbm=cell(d,1);
mat=y{d};
ind_right=cell(d,1);
[mat,dmp]=qr(mat,0);
[ind]=maxvol2(mat);
ind_right{d-1}=ind;
r1=mat(ind,:);
ind_small=cell(d,1);
y{d}=mat/r1;
ind_small{d}=ind;
sbm{d}=r1;
for i=(d-1):-1:2
    core=y{i};
    ncur=size(core,1);
    r2=size(core,2);
    r3=size(core,3);
    core=ten_conv(core,3,r1');
    core=permute(core,[3,1,2]);
    core=reshape(core,[r3*ncur,r2]);
    [core,dmp]=qr(core,0);
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
    ind_small{i}=ind;
    sbm{i}=r1;
    y{i}=core/r1;
    rnew=min(ncur*r3,r2);
    y{i}=reshape(y{i},[r3,ncur,rnew]);
    y{i}=permute(y{i},[2,3,1]);
end
yold=[];
not_converged = true;

while ( swp < nswp && not_converged )
    
    
    
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
    %Compute its SVD, new rank estimate, and new first core
    b=reshape(b,[n1,n2*ra2]);
    
   
    [u,s,v]=svd(b,'econ');
    s=diag(s); nrm=norm(s);
    r=my_chop2(s,eps*nrm);
    r=max([rmin,r]);
    if ( verb ) 
        fprintf('We can push rank %d to %d \n',1,r);  
    end
    u=u(:,1:r);    
    ind=maxvol2(u);
    ind_left{1}=ind;
    u=u/u(ind,:);
    y{1}=u;
%    ind_small{1}=ind;
%     u=u(:,1:r)*diag(s(1:r));
%     v=v(:,1:r);
%     v=reshape(v,[n2,ra2,r]);
%     v=permute(v,[1,3,2]);
%     y{1}=u;
%     y{2}=v;
%norm(vv(:)-full_arr(:))
%Test interpolation property
%keyboard;

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
      r=my_chop2(s,eps*norm(s));
      r=max([r,rmin]);
      %  r=min(r,size(s,1));
  %if ( 
  if ( verb )
    fprintf('We can push rank %d to %d \n',i,r);
  end
  u=u(:,1:r);
  u=reshape(u,[ra1*n1,r]);
  rnew=min(ra1*n1,r);
  [u,rm]=qr(u,0);
  ind=maxvol2(u);
  ind_small{i}=ind;
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
    r=my_chop2(s,eps*norm(s));
    r=max([r,rmin]);
    %r=min(r,size(s,1));
    if ( verb )
      fprintf('We can push rank %d to %d \n',d-1,r);
    end
    u=u(:,1:r);
    v=v(:,1:r);
    s=s(1:r);
    u=u*diag(s);
    
    [ind]=maxvol2(v);
    ind_small{d}=ind;
    mt=v(ind,:);
    u=u*mt';
    u=reshape(u,[ra1,n1,r]);
    u=permute(u,[2,1,3]);
   
    y{d-1}=u;
    y{d}=v/mt;
    ind_right{d-1}=ind;
    v=v/v(ind,:);
    y{d}=v;
    %Here we start iteration backwards
    %this means, that previous step is not that necessary
  for i=d-2:-1:2
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
       %Now we have to fix SVD section 
       [u,s1,v]=svd(b,'econ'); s1=diag(s1);
       r=my_chop2(s1,eps*norm(s1));
       r=max([r,rmin]); 
       if ( verb )
         fprintf('We can push rank %d to %d \n',i,r);
       end
       v=v(:,1:r);
       
  v=reshape(v,[n2*ra2,r]);
  rnew=min(n2*ra2,r);
  ind=maxvol2(v); 
  ind_small{i+1}=ind;
  ind_old=ind_right{i+1};
  ind_new=zeros(d-i,rnew);
  %fprintf('i=%d \n'); 
  for s=1:rnew
     f_in=ind(s);
     w1=tt_ind2sub([n2,ra2],f_in);
     rs=w1(2); js=w1(1);
     ind_new(:,s)=[js,ind_old(:,rs)'];
  end
    ind_right{i}=ind_new;
    mt=v(ind,:);
    s1=s1(1:r);
    u=u(:,1:r)*diag(s1)*mt';
    v=v/mt;
    v=reshape(v,[n2,ra2,rnew]); v=permute(v,[1,3,2]);
    y{i+1}=v;
    u=reshape(u,[ra1,n1,rnew]); u=permute(u,[2,1,3]);
    y{i}=u;
 end
%Fix the first core (separatedly)
%           vv=tt_to_full(y);
% vv=vv(:);
% %Test interpolation property
% k=1;
% for k=1:d-2
%   for s1=1:size(ind_left{k},2);
%      for s2=1:size(ind_right{k+1},2)  
% p0=[ind_left{k}(:,s1)',1,ind_right{k+1}(:,s2)'];
% vv=vv(:);
% p1=tt_sub2ind(sz,p0);
% full_arr=full_arr(:);
% full_arr(p1)-vv(p1)
%      end
%   end
% end
% keyboard;
  
if (isempty(yold) )
  yold=y;
else
 %er=tt_dist2(yold,y)/sqrt(tt_dot(yold,yold)*tt_dot(y,y));
 er=norm(tt_tensor(yold)-tt_tensor(y))/norm(tt_tensor(y));
 yold=y;
 fprintf('er=%3.2e \n',er);
end
   swp=swp+1;
end
% %First core: we compute B(i1,i2,ra)
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
    %Compute its SVD, new rank estimate, and new first core
    b=reshape(b,[n1,n2*ra2]);
    
    [u,s,v]=svd(b,'econ');
    s=diag(s); nrm=norm(s);
    r=my_chop2(s,eps*nrm);
    r=max([rmin,r]);
    %r=min(r,size(s,1));
    if ( verb ) 
        fprintf('We can push rank %d to %d \n',1,r);  
    end
    u=u(:,1:r)*diag(s(1:r));
    v=v(:,1:r);
    v=reshape(v,[n2,ra2,r]);
    v=permute(v,[1,3,2]);
    y{1}=u;
    y{2}=v;
 
%vv=tt_to_full(y);
%norm(vv(:)-full_arr(:))

%keyboard;
return
end
