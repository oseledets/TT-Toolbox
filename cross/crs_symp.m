% This subroutine computes QTT-decomposition of a 
% Characteristic function of a K-dimensional symplex discretized 
% on uniform grid with 2^d points in each direction
K=2;
d=5;
a=0;
b=1;
n=2^d; h=(b-a)/(n-1);
prm=1:K*d; prm=reshape(prm,[d,K]); prm=prm'; prm=reshape(prm,[1,d*K]);
iprm=1:K*d; iprm(prm)=1:K*d;
%Special full tensor parameters
d1=d*K;
tt_old=cell(d1,1);
sz=2*ones(1,d1);
%Describe the array parameters
arr=cell(8,1);
arr_tmp=tt_random(sz,d1,2);
%arr_tmp=tt_zeros(d,sz);
%nrm=sqrt(tt_dot(arr_tmp,arr_tmp));
%arr_tmp=tt_scal(arr_tmp,1./nrm);
arr{1} = d1;
arr{2} = sz;
arr{3} = K;
arr{4} = d; 
arr{5} = h;
arr{6} = a;
arr{7} = iprm;
arr{8}= arr_tmp;

p=tt_crs(arr,1e-12,10);
%p1=tt_crs4(arr,1e-9,20);

