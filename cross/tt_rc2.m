function [y]=tt_rc2(d,n,elem_fun,eps,varargin)
%[Y]=TT_RC2(D,N,ARR,ELEM_FUN,EPS,[OPTIONS])
%The TT-Renormalization-Cross algorithm
%input is: the dimension of the array d, the 
%size vector n, the function that computes the prescribed element
%of an array, accuracy parameter eps
%and also additional arguments specified as 'rmax',10,'nswp',10
%'verb',true and so on
% 'x0', x0
%The elem function has syntax elem_fun(ind); all other additional
%information
%has to be passed as "anonymous" arguments

%The algorithm. It needs left & right index sets to be initialized.
%It needs supecores to be precomputed
%It is a good idea to store also the supercores (there will be d-1,
%i think) & update them; also certain indices maybe removed from the
%index set?? we will not do that; 
%Two possibilities: current index contains maxvol submatrix;
%If not, increase, increase & increase (maybe even beyond the rank!)
%rank one subtraction can be made very stable


%Default parameters
if ( numel(n) == 1 )
  n=n*ones(1,d);
end
nswp=10;
rmax=1000;
x0=[];
verb=true;
vec=false; %If there is a vectorization in options, i.e. if elem_fun
%supports vectorized computations of many elements at once
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
ry=ones(d+1,1); %Initial ranks
%Find fiber maximum
ind=2*ones(d,1);
ind=find_fiber_maximum(elem_fun,ind);
for i=1:d
  i_left{i}=ind(i);
  i_right{i}=ind(i);
  r_left{i}=1;
  r_right{i}=1;
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

%Through all the s-u-p-e-r-c-o-r-e-s
swp=1;
dir='lr';
i=1;
while ( swp <= nswp )
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
  % if ( numel(i1)==2)
  %     if ( i1(1) == 2 && i1(2) == 2 )
  %   keyboard;
  %     end
  % end
  [i1,i2]=new_cross(cur_core,i1,i2,eps); 
       %Increase i_left{i} & i_right{i+1}
       [radd_left,iadd_left]=ind2sub([ry(i);n(i)],i1);
       [iadd_right,radd_right]=ind2sub([n(i+1);ry(i+2)],i2);
       i_left{i}=iadd_left;
       r_left{i}=radd_left;
       %if ( i == 8 )
       %    keyboard
       %end
       i_right{i+1}=iadd_right;
       r_right{i+1}=radd_right;
       ry(i+1)=numel(i1);
       %Zaglushka
       ind_left=get_multi_left(i_left,r_left,ry);
       
       ind_right=get_multi_right(i_right,r_right,ry);
       
      %   fprintf('Index set in question: \n');
      %   ind_left{3}
      % if ( i == 3 && swp == 2 )
      %   keyboard;
      % end
       %keyboard;
       if ( i > 1 ) 
         %Recompute supercore{i-1}  
         super_core{i-1}=compute_supercore(i-1,elem_fun,d,n,ry,ind_left,ind_right,vec);
       end
       if ( i < d - 1)
         super_core{i+1}=compute_supercore(i+1,elem_fun,d,n,ry,ind_left,ind_right,vec);
         %Recompute supercore{i+1}; 
       end

   %Convert iadd1 to i_left and i_right format; compute new
   %ind_left{i+1},ind_right{i+1}
   if ( strcmp(dir,'lr') )
     if ( i == d-1 )
        dir='rl';
     else
        i=i+1;
     end
   else
     if ( i == 1 )
        dir ='lr';
        swp=swp+1;
     else
        i=i-1;
     end
   end
   fprintf('Step %d sweep %d rank=%3.1f \n',i,swp,mean(ry));
end
fprintf('Done! \n');
keyboard;







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

function [i1,i2]=new_cross(mat,i1,i2,eps)
%[i1,i2]=enlarge_cross(mat,i1,i2,eps) 
%Tests whether the current index set gives a reasonable approximation
%to the matrix, and enlarges the basis if necessary.
%The method acts as follows. 

%A-C*A11^{-1} R

%sbm=mat(i1,i2); 
%[u,s,v]=svd(sbm,'econ');
%nrm=norm(mat);s=diag(s); 
%r=numel(find(
%mat_save=mat;
i1=[];
i2=[];
desirata=eps*norm(mat,'fro');
er=norm(mat,'fro');
while ( er > desirata )
   %Find maximal element in mat
   [~,ind]=max(abs(mat(:)));
   ind=tt_ind2sub(size(mat),ind);
   i=ind(1); j=ind(2);
    i1=[i1;i];
    i2=[i2;j];
    u1=mat(:,j); u2=mat(i,:);
    u1=u1/mat(i,j);
    mat=mat-u1*u2;
    er=norm(mat,'fro');
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
  for i=2:d
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
  for i=d-1:-1:1
      %ind_cur=zeros(
      %ind_cur is an array of size ry(i) x i
      ind_cur=zeros(ry(i),d-i+1);
      r_prev=r_right{i};
      ind_prev=ind_right{i+1};
      i_prev=i_right{i};
      for s=1:ry(i)
         ind_cur(s,:)=[i_prev(s), ind_prev(r_prev(s),:)];
      end
      ind_right{i}=ind_cur;
  end

return
end

function [ind]=find_fiber_maximum(elem_fun,ind)
%[ind]=find_fiber_maximum(elem_fun,ind)
%Simple ALS method to compute (approximate) maximum in a tensor
%In fact, need some non-zero
%fact=2; %If the new one <= fact times larger than the previous, stop

return
end
