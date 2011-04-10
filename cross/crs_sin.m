% This subroutine computes QTT-decomposition of a 
% Characteristic function of a K-dimensional symplex discretized 
% on uniform grid with 2^d points in each direction

d=77;
a=0;
b=10^10r;
n=2^d; h=(b-a)/(n-1);
d1=d;
tt_old=cell(d1,1);
sz=2*ones(1,d1);
%Describe the array parameters
arr=cell(7,1);
arr{1} = d1;
arr{2} = sz;
arr{4} = d; 
arr{5} = h;
arr{6} = a;

%p=tt_crs2(arr,1e-9,5);
p1=tt_crs(arr,1e-8,20);
tmp=tt_ones(d,2);
val=tt_conv(p1,tmp);
%val*h

