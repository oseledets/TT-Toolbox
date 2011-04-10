%Prepare a simple test for A*x compression
% A = (I+I), x is vector of all ones
clear all;
d=32;
n=2^d;
%x=tt_ones(2,d);
%vec=1:n; x=vec_to_nd(vec,1e-8,d);
m=4;
%sz=(1:d)*m;
sz=2;
rk=30;
x=tt_random(sz,d,rk);
%a=tt_add(tt_eye(sz,d),tt_eye(sz,d));
a=tt_random(sz,d,1);
for i=1:5
  a = tt_add(a,a);
end
a=tt_vec_to_diag(a);
fprintf('rank of a: %d rank of x: %d \n',tt_erank(a),tt_erank(x));
tic; mvk=tt_mvk(a,x,1e-5); tmvk1=toc;
tic; p=tt_mv(a,x); t_mv=toc; tic; p=tt_compr2(p,1e-9); t_svd=toc;
tic; pp=tt_mvals(a,x,p,5); t_als_5=toc;
tic; pp1=tt_mvals(a,x,tt_random(tt_size(x),d,40),5); t_als_40=toc;
tic; pp2=tt_mvk2(a,x,1e-6); t_newk=toc;
ermvk1=tt_dist2(mvk,p)/sqrt(tt_dot(p,p));
er_als_5=tt_dist2(pp,p)/sqrt(tt_dot(p,p));
er_als_40=tt_dist2(pp1,p)/sqrt(tt_dot(p,p));
er_mvk2=tt_dist2(pp2,p)/sqrt(tt_dot(p,p));
fprintf('SVD:    \ttime=%3.1f\ter=%3.2e\t\nKrylov1:\ttime=%3.1f\ter=%3.2e\nALS_P:     \ttime=%3.1f\ter=%3.2e\t\nALS_R40:\ttime=%3.1f\ter=%3.2e\t\nKrylov2:\ttime=%3.1f\ter=%3.2e\n',t_svd,1e-9,tmvk1,ermvk1,t_als_5,er_als_5,t_als_40,er_als_40,t_newk,er_mvk2);
fprintf('Tschussi! \n');
