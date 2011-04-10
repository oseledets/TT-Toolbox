rmax=40;
eps=1e-6;
n0=2;
nrm_old=0;
d1=3; n=2^d1;
a=0; b=20;
h=(b-a)*1.0/(n-1); 
x=a+(0:n-1)*h;
arr=cell(10,1);

arr{1}=d;
sz=n*ones(1,d);
arr{2}=sz;
arr{3}=x;

p=tt_crs(arr,1e-6,10);
p1=tt_crs(arr,1e-6,20);
tt_dist2(p,p1)/sqrt(tt_dot(p,p))