function [y]=tt_rc(d,n,elem_fun,eps,varargin)
%[Y]=TT_RC(D,N,ARR,ELEM_FUN,EPS,[OPTIONS])
%The TT-Renormalization-Cross algorithm
%input is: the dimension of the array d, the 
%size vector n, the function that computes the prescribed element
%elem_fum
%Default parameters
%nswp=10 (number of DMRG sweeps )
%rmax=1000 (maximal rank of the solution)
%x0=[];  (initial approximation)
%verb = [0,1,2] (verbosity level)
%vec=false (vectorization of the element function)

if ( numel(n) == 1 )
  n=n*ones(d,1);
end
nswp=50;
rmax=1000;
x0=[];
verb=2;
vec=false; %If there is a vectorization in options, i.e. if elem_fun
%supports vectorized computations of many elements at once
change_dir_on=false;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
            rmax=lower(varargin{i+1});
        case 'x0'
            x0=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'vec'
            vec=varargin{i+1};
        case 'change_dir_on'
            change_dir_on=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

%The initial setup is done. 

%We need initial index sets. Say, the vector of all ones (he-he)
%i_left;
%r_left;
%i_right;
%r_right;
i_left=cell(d,1);
i_right=cell(d,1);
r_left=cell(d,1);
r_right=cell(d,1);

if ( isempty(x0) )
ry=ones(d+1,1); %Initial ranks
%Find fiber maximum
ind=2*ones(d,1);
%ind=find_fiber_maximum(elem_fun,ind);
%ind=randi([1,2],[d,1]);
ind=find_fiber_maximum(elem_fun,n,ind);
for i=1:d
  i_left{i}=ind(i);
  i_right{i}=ind(i);
  r_left{i}=1;
  r_right{i}=1;
end
else
  ry=x0.r;
  ps=x0.ps;
  cr=x0.core;
  rm=1;
  x0_old=x0;
  for i=1:d-1
     cur_core=cr(ps(i):ps(i+1)-1);
     cur_core=reshape(cur_core,[ry(i),n(i)*ry(i+1)]);
     cur_core=rm*cur_core;
     cur_core=reshape(cur_core,[ry(i)*n(i),ry(i+1)]);
     [cur_core,rm]=qr(cur_core,0);
     ind=maxvol2(cur_core); ind=ind';
     sbm=cur_core(ind,:);
     cur_core=cur_core / sbm; 
     %q*inv(sbm)*sbm*rm
     rm=sbm*rm;
     cr(ps(i):ps(i+1)-1)=cur_core(:);
     %Convert ind to i_left r_left format
     [radd,iadd]=ind2sub([ry(i);n(i)],ind);
     i_left{i}=iadd;
     r_left{i}=radd;
  end
  cur_core=cr(ps(d):ps(d+1)-1);
  cur_core=reshape(cur_core,[ry(d),n(d)*ry(d+1)]);
  cur_core=rm*cur_core;
  cr(ps(d):ps(d+1)-1)=cur_core(:);
  %And backwards
  rm=1;
  for i=d:-1:2
     cur_core=cr(ps(i):ps(i+1)-1);
     cur_core=reshape(cur_core,[ry(i)*n(i),ry(i+1)]);
     cur_core=cur_core*rm;
     cur_core=reshape(cur_core,[ry(i),n(i)*ry(i+1)]);
     cur_core=cur_core.';
     [cur_core,rm]=qr(cur_core,0);
     ind=maxvol2(cur_core);  ind=ind';
     sbm=cur_core(ind,:);
     cur_core=cur_core/sbm;
     rm=sbm*rm;
     cur_core=cur_core.';
     cr(ps(i):ps(i+1)-1)=cur_core;
     rm=rm.';
     [iadd,radd]=ind2sub([n(i);ry(i+1)],ind);
     i_right{i}=iadd;
     r_right{i}=radd;
  end
  cur_core=cr(ps(1):ps(2)-1);
  cur_core=reshape(cur_core,[ry(1)*n(1),ry(2)]);
  cur_core=cur_core*rm;
  cr(ps(1):ps(2)-1)=cur_core(:);
  x0.core=cr(:);
  %keyboard;
end

%Before we get the multiindex sets (helps to compute s-u-p-e-r-c-o-r-e-s)

ind_left=get_multi_left(i_left,r_left,ry);
ind_right=get_multi_right(i_right,r_right,ry);

%Now we have to compute s-u-p-e-r-c-o-r-e-s (not cores!) using ind_left &
%ind_right
super_core=cell(d-1,1);
for i=1:d-1 %There are d-1 s-u-p-e-r-c-o-r-e-s
   %Compute the index set for the i-th one
   index_set=zeros(d,ry(i),n(i),n(i+1),ry(i+2));
   if ( i == 1 )
       ileft=zeros(1,0);
   else
      ileft=ind_left{i-1};    
   end
   if ( i == d - 1 )
      iright = zeros(1,0);
   else
      iright=ind_right{i+2};
   end
   for s1=1:ry(i)
      for s2=1:ry(i+2)
        for i1=1:n(i)
           for i2=1:n(i+1)
              index_set(:,s1,i1,i2,s2)=[ileft(s1,:),i1,i2,iright(s2,:)];
           end
        end
      end
   end
   M=ry(i)*n(i)*n(i+1)*ry(i+2);
   index_set=reshape(index_set,[d,M]);
   if ( vec ) 
      super_core{i}=reshape(elem_fun(index_set),[ry(i),n(i),n(i+1),ry(i+2)]);
   else
       cur_core=zeros(M,1);
      for k=1:M
         cur_core(k)=elem_fun(index_set(:,k));
      end
      super_core{i}=reshape(cur_core,[ry(i),n(i),n(i+1),ry(i+2)]);
   end
   
end
%Initialization is done. 
%Now to the main iteration
restart=true;
while ( restart )
    
%Through all the s-u-p-e-r-c-o-r-e-s
swp=1;
dir='lr';
i=1;
converged=false;
rank_increase=false;
while ( swp <= nswp && ~converged)
   %The method acts as follows.
   %1) Test if the current supercore is well approximated by the
   %current cross
   %2) If yes, do nothing. 
   %3) If no, add the required indices to the cross and perform
   %modifications of the s-u-p-e-r-c-o-r-e-s
   
   %Nikrena ne soobrajau
   %Bilo: (i1,i2,i3,i4,i5) 
   %Pofiksili (i1,i2), uvelichili r2. Rashirilos' mnojestvo (i2,i3,i4,i5)
   %(nikomu ne nujno)
   %V seredine: (i1,i2) (i3,i4,i5)+ - uvelichilos' mojectvo
   %Kogda uvelichili (i3,i4,i5) viros ry(3), uvelichilos 1 superyadro
   cur_core=super_core{i}; 
   %we have to convert the pair i_left{i} & r_left{i} to the 
   %index set for the matrix cur_core
   
   cur_core=reshape(cur_core,[ry(i)*n(i),n(i+1)*ry(i+2)]);
   %i1=get_full_ind(i_left{i},r_left{i});
   %i2=get_full_ind(i_right{i+1},r_right{i+1});
    %ry(i)*n(i)
    
    i1=r_left{i}+(i_left{i}-1)*ry(i); 
    i2=i_right{i+1}+(r_right{i+1}-1)*n(i+1);
 
  %end
  %if ( swp == 3 && i == 9 )
  %  keyboard;
  %end
      %If the cross is significantly "bad" then we can throw out certain
      %indices and redo, either recompute the full cross
      [iadd1,iadd2]=enlarge_cross(cur_core,i1,i2,eps); 
      if ( ~isempty(iadd1) ) %There will be new elements in supercore{i+1} and supercore{i-1}!
       rank_increase=true;
       %Increase i_left{i} & i_right{i+1}
       [radd_left,iadd_left]=ind2sub([ry(i);n(i)],iadd1);
       [iadd_right,radd_right]=ind2sub([n(i+1);ry(i+2)],iadd2);
       i_left{i}=[i_left{i};iadd_left];
       r_left{i}=[r_left{i};radd_left];
       %if ( i == 8 )
       %    keyboard
       %end
       i_right{i+1}=[i_right{i+1};iadd_right];
       r_right{i+1}=[r_right{i+1};radd_right];
       rold=ry(i+1);
       ry(i+1)=ry(i+1)+numel(iadd1);
       %Zaglushka
       ind_left=get_multi_left(i_left,r_left,ry);
       
       ind_right=get_multi_right(i_right,r_right,ry);
       
       if ( i > 1 ) 
         %Recompute supercore{i-1};
         %supercore{i-1}=zeros(ry(i-1),n(i-1),n(i),ry(i+1))
         radd=numel(iadd1);
         rtmp=ry;  rtmp(i+1)=radd;
         %tmp1=ind_left;
         %tmp1{i}(:,end-radd+1:end)=[];
         tmp2=ind_right;
         tmp2{i+1}(1:rold,:)=[];
         core_add=compute_supercore(i-1,elem_fun,d,n,rtmp,ind_left,tmp2,vec);
         new_core=zeros(ry(i-1),n(i-1),n(i),ry(i+1));
           new_core(:,:,:,1:rold)=super_core{i-1};
         new_core(:,:,:,rold+1:end)=core_add;
         super_core{i-1}=new_core;
         %p1=compute_supercore(i-1,elem_fun,d,n,ry,ind_left,ind_right,vec);
         %keyboard;
       end
       if ( i < d - 1)
         radd=numel(iadd1);
         rtmp=ry;  rtmp(i+1)=radd;
         tmp1=ind_left;
         tmp1{i}(1:rold,:)=[];
         core_add=compute_supercore(i+1,elem_fun,d,n,rtmp,tmp1,ind_right,vec);
         new_core=zeros(ry(i+1),n(i),n(i+1),ry(i+3));
         if ( rold > 0 ) 
           new_core(1:rold,:,:,:)=super_core{i+1};
         end
         new_core(rold+1:end,:,:,:)=core_add;
         super_core{i+1}=new_core;
         %Recompute supercore{i+1}; 
       end
       
       %Now one can check, if which supercore is better approximated
   
       
   end
       if ( i ~= 1 && i ~= d - 1 && change_dir_on) 
         sm=super_core{i-1}; sp=super_core{i+1};
          i1m=r_left{i-1}+(i_left{i-1}-1)*ry(i-1); 
          i2m=i_right{i}+(r_right{i}-1)*n(i);
           i1p=r_left{i+1}+(i_left{i+1}-1)*ry(i+1); 
           i2p=i_right{i+2}+(r_right{i+2}-1)*n(i+2);
           
           sm=reshape(sm,[ry(i-1)*n(i-1),n(i)*ry(i+1)]);
           sp=reshape(sp,[ry(i+1)*n(i+1),n(i+2)*ry(i+3)]);
           smb=sm(:,i2m); [smb,dmp]=qr(smb,0);
           sub_m=smb(i1m,:);
           sm_appr=(smb / sub_m) * sm(i1m,:);
           
           spb=sp(:,i2p); [spb,dmp]=qr(spb,0);
           sub_p=spb(i1p,:);
           sp_appr=(spb / sub_p) * sp(i1p,:);


       end
   %Convert iadd1 to i_left and i_right format; compute new
   %ind_left{i+1},ind_right{i+1}
   
   if ( verb >= 1 ) 
         fprintf('Step %d sweep %d ry(%d): %d -> %d rank=%3.1f \n',i,swp,i+1,ry(i+1)-numel(iadd1),ry(i+1),mean(ry));
   end
   if ( i ~= 1 && i ~= d-1 && change_dir_on )
      er1=norm(sm-sm_appr,'fro')/norm(sm,'fro');
      er2=norm(sp-sp_appr,'fro')/norm(sp,'fro');
   else
     er1=0; er2=0;
   end
       

   if ( (strcmp(dir,'lr') && er1 > er2) && er1 > eps && change_dir_on )
        dir='rl';
        if ( verb >= 1 ) 
                fprintf('i=%d CHANGED DIRECTION LR->RL er1=%3.2e, er2=%3.2e\n',i,er1,er2);
        end

        i=i-1;
        rank_increase=true;
   elseif ( strcmp(dir,'lr')  )
        i=i+1;
   elseif ( (strcmp(dir,'rl') && er2 > er1 ) && er2 > eps && change_dir_on)
        dir='lr';
        i=i+1;
        if ( verb >= 1 ) 
          fprintf('i=%d CHANGED DIRECTION RL->LR er1=%3.2e, er2=%3.2e\n',i,er1,er2);
        end
        rank_increase=true;
   elseif ( strcmp(dir,'rl') )
        i=i-1;
   end
               
           %fprintf('er1=%3.2e, er2=%3.2e \n',er1,er2);
   
   if ( strcmp(dir,'lr') )
     if ( i == d-1 )
        dir='rl';
     end
   else
     if ( i == 1 )
        dir ='lr';
        swp=swp+1;
        if ( rank_increase ) 
          rank_increase=false;
        else
          converged=true;
        end
     end
   end
end
if ( verb >= 1 ) 
fprintf('Converged in %d sweeps  \n',swp);
end
%Now compute the approximation itself. Will it be easy?
%Simple idea: compute cores out of the supercores; they-are-the-same


%% Computation of the approximation 
%Supercore->submatrix->SVD (?) -> new core
ps=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
cr=zeros(ps(d+1)-1,1);

for i=1:d-1
  score=super_core{i};
  score=reshape(score,[ry(i)*n(i),n(i+1)*ry(i+2)]);
  i2=i_right{i+1}+(r_right{i+1}-1)*n(i+1);
  i1=r_left{i}+(i_left{i}-1)*ry(i); 

  %Now do just ordinary check
  r=0; 
  er=2*eps;
  %sbm=score(i1,i2);
  %score=score(:,i2); 
    %We have to solve (approximattttely) the equation
   %Newcore*sbm=score
   %sbm can be rank-deficient
   
  %Try and test the rank
  sbm=score(i1,i2);
  [u,s,v]=svd(sbm,'econ'); s=diag(s);
  %
  maxr=max(size(sbm));
  while ( r < maxr && er > eps )
     r=r+1;
      u0=u(:,1:r); s0=s(1:r); v0=v(:,1:r);
     iu=maxvol2(u0); 
     iv=maxvol2(v0);
     iu1=i1(iu); iv1=i2(iv); 
     w1=score(:,iv1); [w1,dmp]=qr(w1,0);
     w1=w1 / w1(iu1,:);
     appr=w1*score(iu1,:);
     er=norm(score-appr,'fro')/norm(score,'fro');
  end
  %Now the most delicate part: correctly form the core
  
  %cr(pos:pos+ry(i)*n(i)*ry(i+1)-1)=w1(:);
  w2=zeros(ry(i)*n(i),ry(i+1));
  w2(:,iu)=w1;
  cr(ps(i):ps(i+1)-1)=w2(:);
  %norm(w2*score(i1,:)-score,'fro')
  %keyboard;
  %pos=pos+ry(i)*n(i)*ry(i+1);
  %And fix the (i+1)-th supercore 
  %iv provides the truncation ry(i+1)->r 
  
  
 
end
%The last one 
score=super_core{d-1};
 score=reshape(score,[ry(i)*n(i),n(i+1)*ry(i+2)]);
 i1=r_left{d-1}+(i_left{d-1}-1)*ry(d-1); 
score=score(i1,:);
score=reshape(score,[ry(d),n(d),ry(d+1)]);
cr(ps(d):ps(d+1)-1)=score(:);

y=tt_tensor;
y.core=cr;
y.ps=ps;
y.r=ry;
y.n=n;
y.d=d;
%% RANDOMIZED CHECK
%Now perform a randomized check of the approximation
tmp_fun = @(ind) elem_fun(ind) - y(ind); 
%s=100;
%ind=zeros(d,s);
%for i=1:d
%  ind(i)=randi([1,n(i)],1,1);
%end
ind=ones(1,d);
for s=1:5
  [ind]=find_fiber_maximum(tmp_fun,n,ind);
end
  if ( tmp_fun(ind) > 10*eps*elem_fun(ind) ) 
    %This is a restart
    restart=true;
    %Enlarge the index set by ind
    %We have to parse the index set
    rprev=1;
    i=1;
    
    for i=1:d-1
       i_left{i}=[i_left{i};ind(i)];
       if ( i > 1 ) 
         r_left{i}=[r_left{i};ry(i)+1];
       else
         r_left{i}=[r_left{i};1];    
       end
    end
    for i=2:d
      i_right{i}=[i_right{i};ind(i)];
      if ( i < d ) 
        r_right{i}=[r_right{i}; ry(i+1)+1];
      else
        r_right{i}=[r_right{i}; 1];
      end
    end
    ry(2:d)=ry(2:d)+1;
    %Bolshie zaglushki
    ind_left=get_multi_left(i_left,r_left,ry);
       
    ind_right=get_multi_right(i_right,r_right,ry);
    for i=1:d-1
       super_core{i}=compute_supercore(i,elem_fun,d,n,ry,ind_left,ind_right,vec);
    end
    %keyboard;
  else
     restart=false;
  end
end %The restart cycle
return
end

function [super_core]=compute_supercore(i,elem_fun,d,n,ry,ind_left,ind_right,vec)
%[super_core]=compute_supercore(i,elem_fun,n,ry,ind_left,ind_right)
   index_set=zeros(d,ry(i),n(i),n(i+1),ry(i+2));
   if ( i == 1 )
       ileft=zeros(1,0);
   else
      ileft=ind_left{i-1};    
   end
   if ( i == d - 1 )
      iright = zeros(1,0);
   else
      iright=ind_right{i+2};
   end
   for s1=1:ry(i)
      for s2=1:ry(i+2)
        for i1=1:n(i)
           for i2=1:n(i+1)
              index_set(:,s1,i1,i2,s2)=[ileft(s1,:),i1,i2,iright(s2,:)];
           end
        end
      end
   end
   M=ry(i)*n(i)*n(i+1)*ry(i+2);
   index_set=reshape(index_set,[d,M]);
   if ( vec ) 
      super_core=reshape(elem_fun(index_set),[ry(i),n(i),n(i+1),ry(i+2)]);
   else
       cur_core=zeros(M,1);
      for k=1:M
         cur_core(k)=elem_fun(index_set(:,k));
      end
      super_core=reshape(cur_core,[ry(i),n(i),n(i+1),ry(i+2)]);
   end
   

return
end

function [iadd1,iadd2]=enlarge_cross(mat,i1,i2,eps)
%[iadd1,iadd2,flag]=enlarge_cross(mat,i1,i2,eps) 
%Flag is raised if the initial cross is "very" bad
%Tests whether the current index set gives a reasonable approximation
%to the matrix, and enlarges the basis if necessary.
%The method acts as follows. 

magic_eps=1e-12; %just some magic


%sbm=mat(i1,i2); 
%[u,s,v]=svd(sbm,'econ');
%nrm=norm(mat);s=diag(s); 
%r=numel(find(
%mat_save=mat;
i1save=i1;
i2save=i2;
sbm=mat(i1,i2); 

%Compute the initial approximation. The computation is based on the 
%LU decomposition of the sbm until the submatrix is well approximated
[u,s,v]=svd(sbm,'econ');
s=diag(s); rm=my_chop2(s,min(magic_eps,eps)*norm(s)); 
u=u(:,1:rm); v=v(:,1:rm);
indu=maxvol2(u); indv=maxvol2(v);
i1=i1(indu);i2=i2(indv); 

cbas=mat(:,i2);
rbas=mat(i1,:); 
%rbas=rbas.';
%[cbas,~]=qr(cbas,0);
%[rbas,~]=qr(rbas,0);
%phi=cbas'*mat*rbas;
%appr=cbas*phi*rbas';
[cbas,~]=qr(cbas,0);
qsbm=cbas(i1,:);
% 
% if ( cond(qsbm) > 1e14 )
%   flag = true;
%   appr=zeros(size(mat));
% else
%   mm=cbas / qsbm;
%   if ( max(abs(mm(:))> magic_factor ))
%     flag=true; appr=zeros(size(mat));
%   else
    appr=cbas / qsbm * rbas;    
%  end
%end


iadd1=[];
iadd2=[];

if ( norm(mat(:)-appr(:),'fro') < eps*norm(mat(:),'fro') )
   rnew=numel(i1);
else
   desirata=eps*norm(mat(:),'fro');
   er=norm(mat(:)-appr(:),'fro');
   mat=mat-appr;
   while ( er > desirata )
      %Find maximal element in mat
      [~,ind]=max(abs(mat(:)));
      ind=tt_ind2sub(size(mat),ind);
      i=ind(1); j=ind(2);
      %if ( ~any(i1==i) ) 
        iadd1=[iadd1;i];
      %end
      %if ( ~any(i2==j) )
        iadd2=[iadd2;j];
      %end
      u1=mat(:,j); u2=mat(i,:);
      u1=u1/mat(i,j);
      mat=mat-u1*u2;
      er=norm(mat(:),'fro');
   end
end
% i0=[i1;iadd1]; j0=[i2;iadd2];
% sbm=mat_save(i0,j0);
% ss=svd(sbm);
% fprintf('%10.10e \n',ss)
% sbm=mat_save(i1,i2);
% fprintf('AND PREVIOUS \n');
% ss=svd(sbm);
% fprintf('%10.10e \n',ss)

return
end

function [ind_left]=get_multi_left(i_left,r_left,ry)
%[ind_left]=get_multi_left(i_left,r_left,ry)
%Computes (all) left multiindex sets for a given compact representation
  d=numel(i_left);
  ind_left=cell(d,1);
  ind_left{1}=i_left{1};
  for i=2:d-1
      %ind_cur=zeros(
      %ind_cur is an array of size ry(i-1) x i
      ind_cur=zeros(ry(i+1),i);
      r_prev=r_left{i};
      ind_prev=ind_left{i-1};
      i_prev=i_left{i};
      for s=1:ry(i+1)
          %ind_prev(p1(s
         ind_cur(s,:)=[ind_prev(r_prev(s),:),i_prev(s)];
      end
      ind_left{i}=ind_cur;
  end
return
end

function [ind_right]=get_multi_right(i_right,r_right,ry)
%[ind_right]=get_multi_right(i_right,r_right,ry)
%Computes (all) right multiindex sets for a given compact representation
  d=numel(i_right);
  ind_right=cell(d,1);
  ind_right{d}=i_right{d};
  for i=d-1:-1:2
      %ind_cur=zeros(
      %ind_cur is an array of size ry(i) x i
      ind_cur=zeros(ry(i),d-i+1);
      r_prev=r_right{i};
      ind_prev=ind_right{i+1};
      i_prev=i_right{i};
      for s=1:ry(i)
         % fprintf('i=%d s=%d \n',i,s);
         ind_cur(s,:)=[i_prev(s), ind_prev(r_prev(s),:)];
      end
      ind_right{i}=ind_cur;
  end
  
  
return
end

function [ind]=find_fiber_maximum(elem_fun,n,ind)
%[ind]=find_fiber_maximum(elem_fun,n,ind)
%Simple ALS method to compute (approximate) maximum in a tensor
%In fact, need some non-zero
%fact=2; %If the new one <= fact times larger than the previous, stop
%Compute the fibers; find maximum along them; 
d=numel(n);
mx=elem_fun(ind); mx=abs(mx);
git=2;
for s=1:git
for k=1:d
  ind_tmp=ind;
  for i=1:n(k)
    ind_tmp(k)=i;
    val=abs(elem_fun(ind_tmp));
    if ( val >= mx ) 
       ind(k)=i;
       mx=val;
       %fprintf('mx=%f \n',mx);
    end
  end
end
end
return
end
