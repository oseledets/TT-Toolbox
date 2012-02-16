function [tt]=tt_exp2(mat,eps,N,rmax)
%Computation of the matrix exponential of a matrix in TT format
%   [TT]=TT_EXP2(MAT,EPS,N,RMAX) This function computes matrix exponential
%   using scaling and squaring method of the TT-matrix TT. EPS is the
%   accuracy parameter, N is the number of summand for the local Taylor
%   series (N=10 is usually enough). RMAX is the TT-rank bound and TT_MAT
%   is the result
%
%
% TT-Toolbox 2.2, 2009-2011
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
if ( nargin == 3) 
  rmax=1000;
end
nrm=abs(norm(mat,2));

n0=floor(max(floor(log2(10*nrm)),0))+1;
w0=mat*(1.0/2^n0);
sz=size(mat); sz=sz(:,1);
e=tt_eye(sz,ndims(mat)); e=tt_matrix(e);
tt=e;

for k=(N-1):-1:1
    %fprintf('Iteration %d/1 \n',k);
    %if ( tt_erank(tt_mat)*tt_erank(w0) > rmax )
    %  tt_mat=ttm_add(tt_scal(tt_mmals(tt_mat,w0,xs,3),1.0/k),e);  
    %else
    %  tt_mat=ttm_add(tt_scal(tt_mm(tt_mat,w0),1.0/k),e);  
    %end
    tt=e+(tt*w0)/k;
    tt=round(tt,eps);
    %tt_mat=tt_mat_compr(tt_mat,eps);
end
%keyboard;
for k=1:n0
   % reff=sqrt(tt_mem(tt_mat)/(n*n*d));
   %fprintf('Iteration %d/%d Rank=%3.2f\n',k,n0,reff);
   tt=round(tt*tt,eps,rmax);
   %else
   %  tt_mat=tt_mm(tt_mat,tt_mat);
   %end
   %tt_mat=tt_mat_compr(tt_mat,eps);
   %fprintf('nrm=%3.2e \n',tt_dot(tt_mat1,tt_mat1));
end

return
end
