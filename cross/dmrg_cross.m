function [y]=dmrg_cross(d,n,fun,eps,varargin)
%DMRG-cross method for the approximation of TT-tensors
%   [A]=DMRG_CROSS(D,N,FUN,EPS,OPTIONS) Computes the approximation of a
%   given tensor via the adaptive DMRG-cross procedure. The input is a pair
%   (D,N) which determines the size of the tensor (N can be either a
%   number, or array of mode sizes). FUN is the function to compute
%   a prescribed element of a tensor (FUN(IND)), or it can be vectorized to
%   compute series of elements of a tensor (see OPTIONS) To pass parameters 
%   to FUN please use anonymous function handles. EPS is the accuracy 
%   of the approximation.Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following) 
%   The list of option names and default values are:
%       o nswp - number of DMRG sweeps [10]
%       o vec  - Fun is vectorized [ true | {false} ]
%       o verb - output debug information [ true | {false} ]
%       o y0   - initial approximation [random rank-2]
%       o radd - minimal rank change [0]
%       o rmin - minimal rank that is allows [1]
%       
%Difference to TT_CRS4: vectorized computation of subtensors added
%TT_CROSS_ELEM is a function handle than computes a set of elements 
%of array A

%Default parameters
rmin=1;
verb=true;
radd=0;
nswp=10;
y=[];
vec=false;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'y0'
            y=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'rmin'
            rmin=varargin{i+1};
        case 'radd'
            radd=varargin{i+1};
        case 'vec'
            vec=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if ( numel(n) == 1 )
   n=n*ones(d,1);
end
yold=[];

sz=n;
if (isempty(y) )
    y=tt_rand(sz,d,2); 
end
%y=round(y,0); %To avoid "overranks"
ry=y.r;

%Here we need: "left indices" and right indices; We consequently 
%orthogonalize y from left to right (assuming that from right-to-left 
%they are already orthogonalized)
%At this point we assume that: Y is orthogonal from right-to-left
%and ind_right correspond to correct max-vol positions in this matrix

%Warmup procedure: orthogonalization from right to left of the initial
%approximation & computation of the index sets
ind_left=cell(d,1);
er=2*eps;
swp=1;

r1=1;
for i=d:-1:2
    cr=y{i}; cr=reshape(cr,[ry(i)*n(i),ry(i+1)]);
    cr = cr*r1; cr=cr.';
    [cr,~]=qr(cr,0);
    [core,~]=qr(core,0);
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
    megaindex=zeros(d,n1,n2,ra2);
    for i=1:n1
      for j=1:n2
         for s2=1:ra2
             ind=[i,j,ind2(:,s2)'];
            megaindex(:,i,j,s2)=ind;
         end
      end
    end
    megaindex=reshape(megaindex,[d,n1*n2*ra2]);
    b=tt_cross_elem(megaindex,arr);
        
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
       %b=zeros(ra1,n1,n2,ra2);
       %Create MEGAINDEX array
       megaindex=zeros(d,ra1,n1,n2,ra2); %These are all n1*n2*ra1*ra2 elements to compute (?)
       
       for i1=1:n1
         for i2=1:n2
           for s1=1:ra1
             for s2=1:ra2
                ind=[ind1(:,s1)',i1,i2,ind2(:,s2)'];
         megaindex(:,s1,i1,i2,s2)=ind;
               %         b(s1,i1,i2,s2)=(ind,arr);
             end
           end
         end
       end
       %Compute b
       megaindex=reshape(megaindex,[d,ra1*n1*n2*ra2]);
       b=tt_cross_elem(megaindex,arr);
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
%     for i1=1:n1
%       for i2=1:n2
%          for s1=1:ra1
%             ind=[ind1(:,s1)',i1,i2];
%             b(s1,i1,i2)=tt_cross_elem(ind,arr);
%          end
%       end
%     end
%     
    
    
       megaindex=zeros(d,ra1,n1,n2); 
           
       for i1=1:n1
         for i2=1:n2
           for s1=1:ra1
            ind=[ind1(:,s1)',i1,i2];
                 megaindex(:,s1,i1,i2)=ind;
           end
         end
       end
      
       %Compute b
       megaindex=reshape(megaindex,[d,ra1*n1*n2]);
       b=tt_cross_elem(megaindex,arr);
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
       
       megaindex=zeros(d,ra1,n1,n2,ra2); %These are all n1*n2*ra1*ra2 elements to compute (?)
       
       for i1=1:n1
         for i2=1:n2
           for s1=1:ra1
             for s2=1:ra2
                ind=[ind1(:,s1)',i1,i2,ind2(:,s2)'];
         megaindex(:,s1,i1,i2,s2)=ind;
               %         b(s1,i1,i2,s2)=my_elem(ind,arr);
             end
           end
         end
       end
       %Compute b
       megaindex=reshape(megaindex,[d,ra1*n1*n2*ra2]);
       b=tt_cross_elem(megaindex,arr);
       b=reshape(b,[ra1*n1,n2*ra2]);
   
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
%     b=zeros(n1,n2,ra2);
%     for i=1:n1
%       for j=1:n2
%          for s2=1:ra2
%              ind=[i,j,ind2(:,s2)'];
%             b(i,j,s2)=my_elem(ind,arr);
%          end
%       end
%     end
    
     megaindex=zeros(d,n1,n2,ra2); 
           
       for i1=1:n1
         for i2=1:n2
           for s2=1:ra2
             ind=[i1,i2,ind2(:,s2)'];
             megaindex(:,i1,i2,s2)=ind;
           end
         end
       end
      
       %Compute b
       megaindex=reshape(megaindex,[d,n1*n2*ra2]);
       b=tt_cross_elem(megaindex,arr);

    
    
    
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
