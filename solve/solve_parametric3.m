function [x,hst1,hst2,sols]=solve_parametric3(mat,par,rhs0,x0,eps,niter)
%[X] = SOLVE_PARAMETRIC3(MAT,PAR,RHS,EPS,NITER)
% MAT, PAR --- define matrix, please look at examples
% RHS --- right-hand side
% EPS --- accuracy parameter
% NITER --- number of iterations

%M & d0 & p are really not needed at all

%Parameters
solver='maxvol';
%solver='als';
%Initialization

if ( isempty(x0) )
   sols_all=[];
   x=rhs0;
else
   x=x0;
   sols_all=chunk(x,1,1); sols_all=full(sols_all);
   sols_all=squeeze(sols_all);
   [sols_all,~]=qr(sols_all,0);
end

rhs=rhs0;
nsolves=0;
num_mat=numel(mat);
n=size(mat{1},1); %This is the size of the first, "large" dimension

for iter=1:niter
<<<<<<< Updated upstream
 
    %We need to check the residue. Sorry.
    x1=chunk(x,1,1);
   rx=x1.r; rx=rx(2); %Block size
    x2=chunk(x,2,x.d);
    x1=full(x1);
    x1=reshape(x1,[numel(x1)/rx,rx]);
    %y=zeros(n,rx,num_mat);
    y=zeros(n,num_mat,rx);
    %Multiply in the physical space
    for i=1:num_mat
       for j=1:rx
          %y(:,i,j)=mat{i}*x1(:,j);
          y(:,i,j)=mat{i}*x1(:,j);
       end
    end
    y=reshape(y,[n,numel(y)/n]);
    y0=[];
    for i=1:size(y,2)
      y0=[y0,tt_tensor(y(:,i))];
    end
    y=y0;
    %Multiply parametric parts
    %parp=diag(par)*x2;
    %parp=mvk2(diag(par),x2,eps);
    %parp=round(parp,eps);
    %res=kron(y,parp); 
    %res=res-rhs;
    %res=round(res,eps);
    res=rhs;
  er=norm(res)/norm(rhs);
 fprintf('it=%d error=%3.2e\n',iter,er);
if ( strcmp(solver,'als') ) 
 x1=chunk(x,2,ndims(x));
 f1=chunk(rhs,2,ndims(rhs));
 %Determine the ranks and the blocking matrix
 rm=rank(par,1);
 rv=rank(x1,1);
 rf=rank(f1,1);
 
 ff=chunk(rhs,1,1); ff=full(ff); ff=reshape(ff,[numel(ff)/rf,rf]);
 blk=zeros(rm,rv);
 for_rhs=zeros(rf,rv);
 for i=1:rm
   for j=1:rv
      blk(i,j)=dot(row(par,i),row(x1,j));
=======
tic; 
 
 intr=x;
 [~,~,~,~,ind_right]=tt_canform1(tt_compr2(core(intr),1e-8));
 %ind_right{2} is used to compute interpolation points 
 r=size(ind_right{1},2);
 ind_points=zeros(M,r);
 for s=1:r
   ind=ind_right{1}(:,s);
   for i=1:M
       ind0=sub_to_ind(ind((i-1)*d0+1:i*d0),2*ones(1,d0));
       ind_points(i,s)=ind0;
>>>>>>> Stashed changes
   end
 end
 for i=1:rf
   for j=1:rv
     for_rhs(i,j)=dot(row(f1,i),row(x1,j));
   end
 end
 %Solve it rv times
 nx=x.n;
 sols_add=zeros(nx(1),rv);
 for s=1:rv
    bm0=mat{1}*blk(1,s);
    for k=2:rm
       bm0=bm0+mat{k}*blk(k,s);
    end
    rhs1=ff(:,1);
    for k=2:rf
       rhs1=rhs1+ff(:,k)*for_rhs(k,s);
    end
    sols_add(:,s)=bm0 \ rhs1;
    nsolves=nsolves+1;
 end

 
elseif ( strcmp(solver,'maxvol') )
  %This section has to be rewritten. We have to compute qr of x1, its
  %maxvol, and corresponding submatrices in A & u. This is more or
  %less done in the dmrg_eigb. It is better not to pollute the
  %code with unneccesary stuff and write a separate function 
  %instead of the old-fashioned tt_canform1 however, now we can just
  %simply kill the beast by "computing all elements". This would
  %be also helpfull for working with hanging matrices
  
  x1=chunk(x,2,ndims(x));
  f1=chunk(rhs,2,ndims(rhs));
  %Determine the ranks and the blocking matrix
  rm=rank(par,1);
  rv=rank(x1,1);
  rf=rank(f1,1);
 
  ff=chunk(rhs,1,1); ff=full(ff); ff=reshape(ff,[numel(ff)/rf,rf]);

  

  [~,~,~,~,ind_right]=tt_canform1(tt_compr2(core(x),1e-8));
  %ind_right{2} is used to compute interpolation points in par & f1
  ii=ind_right{1};  
  r=size(ii,2);
    nx=x.n;
  sols_add=zeros(nx(1),r);

  for s=1:r
    ind=ii(:,s);
    %Compute the reduction coefficient
    cm=par(ind);
    cf=f1(ind);
    bm0=cm(1)*mat{1};
    for i=2:rm
       bm0=bm0+cm(i)*mat{i};
    end
    rhs1=cf(1)*ff(:,1);
    for i=2:rf
        rhs1=rhs1+cf(i)*ff(:,i);
    end
    sols_add(:,s)=bm0 \ rhs1;
    nsolves=nsolves +  1;

  end
 end
  if ( isempty(sols_all) )
   [sols_all,~]=qr(sols_add,0);
 else
     sols_all=reort(sols_all,sols_add);
 end
 
 %Compute projection matrices
 %mat_small=zeros(r,r,size(mat,2));
 mat_small=[];
 
 r=size(sols_all,2);
 for i=1:num_mat
 %  %mat_small(:,:,i)=sols_all'*(mat{i})*sols_all;
   wm=sols_all'*mat{i}*sols_all; wm=tt_tensor(wm(:));
   mat_small=[mat_small,wm];
 end
 mat_small=tt_matrix(mat_small,r,r);
 rhs_small=ttm(rhs,1,sols_all);
 
 sol_prev=ttm(x,1,sols_all);
 
%Simple dmrg_solve2 solver for the stochastic part
 nm=kron(mat_small,diag(par));
 sol_red=dmrg_parb(mat_small,par,rhs_small,eps/10,sol_prev,[],3,false);
 %keyboard;
 %return
 %sol_red=dmrg_ssolve2(nm,rhs_small,sol_prev,eps,eps,[],3,[],false);
 x=ttm(sol_red,1,sols_all');
 x=round(x,eps);
 tf=toc;    
%We need to check the residue. Sorry.
    x1=chunk(x,1,1);
   rx=x1.r; rx=rx(2); %Block size
    x2=chunk(x,2,x.d);
    x1=full(x1);
    x1=reshape(x1,[numel(x1)/rx,rx]);
    %y=zeros(n,rx,num_mat);
    y=zeros(n,num_mat,rx);
    %Multiply in the physical space
    for i=1:num_mat
       for j=1:rx
          y(:,i,j)=mat{i}*x1(:,j);
       end
    end
    y=reshape(y,[n,numel(y)/n]);
    y0=[];
    for i=1:size(y,2)
      y0=[y0,tt_tensor(y(:,i))];
    end
    y=y0;
    %Multiply parametric parts
    parp=diag(par)*x2;
    parp=round(parp,1e-8);
    res=kron(y,parp); 
    res=res-rhs;
    res=round(res,1e-8);
 
  er=norm(res)/norm(rhs);
 fprintf('it=%d error=%3.2e time=%3.2f rank=%3.1f \n',iter,er,tf,mean(rank(x)));
 
end
fprintf('Total number of solves: %d Accuracy: %3.2e rank=%3.1f \n',nsolves,er);
return
end
