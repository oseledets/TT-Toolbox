function [x,hst1,hst2,sols]=solve_parametric3(mat,M,d0,p,par,rhs0,x0,eps,niter)
%[X] = SOLVE_PARAMETRIC3(MAT,PAR,RHS,EPS,NITER)
% MAT, PAR --- define matrix, please look at examples
% RHS --- right-hand side
% EPS --- accuracy parameter
% NITER --- number of iterations

%Simple iterative solver
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
    parp=diag(par)*x2;
    parp=round(parp,1e-8);
    res=kron(y,parp); 
    res=res-rhs;
    res=round(res,1e-8);
 
  er=norm(res)/norm(rhs);
 fprintf('it=%d error=%3.2e\n',iter,er);
 
 
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
   end
 end
 nx=x.n;
 sols_add=zeros(nx(1),r);
 for s=1:r
   bm0=mat{1};
   for k=2:M+1
     bm0=bm0+p(ind_points(k-1,s))*mat{k};
   end
   %Form right-hand sides for the equation at multiparameter point
   rhs1=tt_elem2(core(rhs0),ind_right{1}(:,s));
   sols_add(:,s)=bm0 \ rhs1;
   nsolves=nsolves +  1;
   %keyboard;
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
 sol_red=dmrg_parb(mat_small,par,rhs_small,eps,sol_prev,[],3,false);
 %keyboard;
 %return
 %sol_red=dmrg_solve2(nm,rhs_small,sol_prev,eps,eps,[],3,[],false);
 x=ttm(sol_red,1,sols_all');
 x=round(x,eps);
end
fprintf('Total number of solves: %d Accuracy: %3.2e \n',nsolves,er);
return
end