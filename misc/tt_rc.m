function [y]=tt_rc(d,n,elem_fun,eps,varargin)
%[Y]=TT_RC(D,N,ARR,ELEM_FUN,EPS,[OPTIONS])
%The TT-Renormalization-Cross algorithm
%input is: the dimension of the array d, the 
%size vector n, the function that computes the prescribed element
%of an array, accuracy parameter eps
%and also additional arguments specified as 'rmax',10,'nswp',10
%'verb',true and so on
% 'x0', x0
%The elem function has syntax elem_fun(ind); all other additional
%information
%has to be pass as "anonymous" arguments

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
for i=1:d
  i_left{i}=1;
  i_right{i}=1;
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

for i=1:d-1
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
    i1=i_left{i}+(r_left{i}-1)*n(i); 
    i2=r_right{i+1}+(i_right{i+1}-1)*ry(i+2);
   [iadd1,iadd2]=enlarge_cross(cur_core,i1,i2,eps); 
   if ( ~isempty(iadd1) ) %There will be new elements in supercore{i+1} and supercore{i-1}!
       %Increase i_left{i} & i_right{i+1}
       [radd_left,iadd_left]=ind2sub([ry(i),n(i)],iadd1);
       [iadd_right,radd_right]=ind2sub([n(i+1),ry(i+2)],iadd2);
       i_left{i}=[i_left{i};iadd_left];
       r_left{i}=[r_left{i};radd_left];
       i_right{i+1}=[i_right{i+1};iadd_right];
       r_right{i+1}=[r_right{i+1};radd_right];
       ry(i+1)=ry(i+1)+numel(iadd1);
       %Zaglushka
       ind_left=get_multi_left(i_left,r_left,ry);
       ind_right=get_multi_right(i_right,r_right,ry);
       %keyboard;
       if ( i > 1 ) 
         %Recompute supercore{i-1}  
         super_core{i-1}=compute_supercore(i-1,elem_fun,d,n,ry,ind_left,ind_right,vec);
       end
       if ( i < d - 1)
         super_core{i+1}=compute_supercore(i+1,elem_fun,d,n,ry,ind_left,ind_right,vec);
         %Recompute supercore{i+1}; 
       end
       keyboard;
   end
   %Convert iadd1 to i_left and i_right format; compute new
   %ind_left{i+1},ind_right{i+1}
end








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
%[iadd1,iadd2,rnew]=enlarge_cross(mat,i1,i2,eps) 
%Tests whether the current index set gives a reasonable approximation
%to the matrix, and enlarges the basis if necessary.
%The method acts as follows. 
%Compute the column basis; compute the row basis; 
cbas=mat(:,i2);
rbas=mat(i1,:); 
%[cbas,~]=qr(cbas,0);
%[rbas,~]=qr(rbas,0);
%phi=cbas'*mat*rbas;
%phi=cbas \ mat / rbas;
%appr=cbas*phi*rbas;
sbm=mat(i1,i2);
appr=cbas / sbm * rbas;
iadd1=[];
iadd2=[];

if ( norm(mat-appr,'fro') < eps*norm(mat,'fro') )
   rnew=numel(i1);
else
   desirata=eps*norm(mat,'fro');
   er=norm(mat-appr,'fro');
   mat=mat-appr;
   n1=size(mat,1); m1=size(mat,2);
   while ( er > desirata )
      %Find maximal element in mat
      [pv,ind]=max(abs(mat(:)));
      ind=tt_ind2sub(size(mat),ind);
      i=ind(1); j=ind(2);
      iadd1=[iadd1,i];
      iadd2=[iadd2,j];
      u1=mat(:,j); u2=mat(i,:);
      u1=u1/mat(i,j);
      mat=mat-u1*u2;
      er=norm(mat,'fro');
   end
end


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
      ind_cur=zeros(ry(i+1),d-i+1);
      r_prev=r_right{i+1};
      ind_prev=ind_right{i+1};
      i_prev=i_right{i+1};
      for s=1:ry(i+1)
          %ind_prev(p1(s
         ind_cur(s,:)=[i_prev(s), ind_prev(r_prev(s),:)];
      end
      ind_right{i}=ind_cur;
  end

return
end
