function [sol1]=als_solve2(mat,rhs,eps,sol,niter)
%[SOL]=ALS_SOLVE2(MAT,RHS,EPS,[SOL],[NITER])
%Computes approximate solution to the linear system MAT*SOL=RHS
%such as the residue is bounded by ||MAT*SOL-RHS||<EPS*||RHS||
%or maximum number of iterations is achieved

if ( nargin < 4 || isempty(sol) )
   sol=zeros(numel(rhs),1);   
end
if ( nargin < 5 || isempty(niter) )
  niter=1000;
end
%Magic parameters
loc_slv=10;
m0=tt_matrix(mat);
nrm=norm(rhs);
rhs1=rhs-m0*sol;
sol1=zeros(size(sol));
i=1;
er=2*eps*nrm;
m=size(mat{2},1);
while ( i < niter && er > eps*nrm )

   [x,y]=als_solve(mat,rhs1,randn(m,1),loc_slv); rhs1=rhs1-(m0*kron(y,x)); 
   sol1=sol1+kron(y,x);
   er=norm(rhs1); 
   fprintf('i=%d, er=%3.2e \n',i,er/nrm);
   i=i+1;
   
end
    

return
end