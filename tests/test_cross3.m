p=4;

L=6; n=2^L-1; m=(n+1)/p; %Step size
a=0;
b=1;
h=(b-a)/(n+1);


%a2=ones(n+1,n+1); %Values of a2 are defined at midpoints
%a2(1:m,1:m)=10^4;
%mat=Fd_mtx2(a2);
%rhs=ones(n*n,1);
%sol=mat \ rhs; 
%mesh(reshape(sol,n,n))
%return
%The solution is defined at points (1:n)/(n+1)
%mat=Fd_mtx2(a2);
%fprintf('after Fd_mtx2 \n');
%keyboard;
%The diffusion coefficient is defined at points (0.5:n+0.5)/(n+1);
xsol=(1:n)*h;
xa=(0.5:n+0.5)*h;

%Generate right coefficient

a2=zeros(n+1,n+1);
mat=[];
for i=1:p
  for j=1:p
      tmp=zeros(n+1,n+1);
      tmp((i-1)*m+1:i*m,(j-1)*m+1:j*m)=1;
      a2=a2+tmp;
      mat{i+(j-1)*p}=Fd_mtx2(tmp); 
  end
end

smin=1;
smax=2;
d=6;
h=(smax-smin)/(2^d-1);
rhs=ones(n*n,1);

fun = @(ind) lars(p,ind,mat,rhs,n,d,smin,h,(n+1)/2,(n+1)/2);
y=tt_rc(p*p*d,2,fun,1e-5,'nswp',40,'change_dir_on',false);
