function [tt_mat]=tt_exp2(tt,eps,N,rmax)
%[TT_MAT]=TT_EXP2(TT,EPS,N,RMAX)
%Computation of the matrix exponential of a matrix in TT format
%Using scaling-and-squaring method
%Input: TT --- matrix in TT format 
%EPS --- accuracy parameter
%N is the number of summands in Taylor series
%(N=10 is usually enough)
%RMAX is the maximal rank for matrix-by-matrix product
%TT_MAT is the result
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
d=size(tt,1);
n=size(tt{1},1); %Warning: only equal sizes are valid here
e=tt_eye(n,d);
xs=tt_random(n*n,d,rmax); 
xs=tt_vec_to_mat(xs,n,n);
nrm=sqrt(tt_dot(tt_mat_to_vec(tt),tt_mat_to_vec(tt)));
n0=floor(max(floor(log2(nrm)),0))+1;
w0=tt_scal(tt,1.0/2^(n0));
tt_mat=e;
for k=(N-1):-1:1
    %fprintf('Iteration %d/1 \n',k);
    if ( tt_erank(tt_mat)*tt_erank(w0) > rmax )
      tt_mat=ttm_add(tt_scal(tt_mmals(tt_mat,w0,xs,3),1.0/k),e);  
    else
      tt_mat=ttm_add(tt_scal(tt_mm(tt_mat,w0),1.0/k),e);  
    end
    tt_mat=tt_mat_compr(tt_mat,eps);
end
%keyboard;
for k=1:n0
   % reff=sqrt(tt_mem(tt_mat)/(n*n*d));
   %fprintf('Iteration %d/%d Rank=%3.2f\n',k,n0,reff);
   if ( tt_erank(tt_mat)^2 > rmax ) 
     tt_mat=tt_mmals(tt_mat,tt_mat,xs,3);
   else
     tt_mat=tt_mm(tt_mat,tt_mat);
   end
   tt_mat=tt_mat_compr(tt_mat,eps);
   %fprintf('nrm=%3.2e \n',tt_dot(tt_mat1,tt_mat1));
end
