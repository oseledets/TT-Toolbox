% This subroutine computes QTT-decomposition of a 
% Characteristic function of a K-dimensional symplex discretized 
% on uniform grid with 2^d points in each direction
%Special full tensor parameters
d1=20;
sz=2*ones(1,d1);
%sz=[2,2,2,2];
%Describe the array parameters
arr=cell(7,1);
arr{1} = d1;
arr{2} = sz;

p=tt_crs(arr,1e-9,5);
p1=tt_crs2(arr,1e-9,10);
p2=tt_crs3(arr,1e-9,10);
