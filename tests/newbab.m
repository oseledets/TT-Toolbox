%% Here the matrix is defined
%% Example from Babushka SIREV paper --- characteristic functions 
%% of 4 circles

M=4;
mat=cell(M+1,1);
n=256; n0=n;
a=0;
b=1;
h=(b-a)/(n+1);
x1=a+(0:n+1)*h;
a1=ones(n+2,n+2);
mat{1}=Fd_mtx(2,a1,0);
a2=zeros(n+2,n+2);
cen=cell(M,1); %Centers
cen{1}=[0.2,0.2];
cen{2}=[0.2,0.8];
cen{3}=[0.8,0.8];
cen{4}=[0.8,0.2];
rad={0.1,0.1,0.1,0.1};
fac={1,0.9,0.75,0.6};
for i=1:M
   tmp=zeros(n+2,n+2);
   for i1=1:n+2
     for j1=1:n+2
        if ( (x1(i1)-cen{i}(1))^2+(x1(j1)-cen{i}(2))^2 <= (rad{i})^2)
            tmp(i1,j1)=1.0;
        else
            tmp(i1,j1)=0.0;
        end
     end
   end
   mat{i+1}=Fd_mtx(2,fac{i}*tmp,0);   
end

mat=mat'; %It should be a ROW MATRIX
%Parameter grid
%d0=8; np=2^d0;
%a0=-0.99; b0=0; h0=(b0-a0)/(np-1); p=a0+(0:np-1)*h0;
d0=2; np=2^d0; % This is more than enough for a ChebGrid
%a0=-1; b0=1; 
a0=-0.99; b0=0;
%a0=0.499; b0=0.5;
%[p,ww1]=lgwt(np,a0,b0);
%Just uniform grid & Simpson rule
hp=(b0-a0)/(np-1);

t=(0:np-1)/(np-1); p=a0+(b0-a0)*t;
ww1=ones(size(p));
%Generate weights for the Simpson rule:
%1 4 2 4 2 4 ... 4 1
%
ww1(1)=hp/2;
ww1(numel(ww1))=hp/2;
sw=true;
pos=2;
while ( pos <= numel(ww1)-1 )
  ww1(pos)=hp;
  pos=pos+1;
end

%e=ones(size(p));
%fprintf('MEAN: %10.10f \n',dot(f,ww1));
%return


%h0=(b0-a0)/(np-1); p=a0+(0:np-1)*h0;

%a0=0; b0=1; h0=(b0-a0)/(np+1); x0=a0+(1:np)*h0;
%x0=norminv(x0);


dm=vec_to_nd(p,1e-8,d0);
e=tt_ones(d0,2);
%Construct full operator
%e=eye(np);
%e=diag(e);

%%

%The approximation scheme is: ``a-la'' Black-box approximation
%for physical dimensions: set of M (sparse) matrices stored in mat
%and M*P*P*P-dimensional tensor of parameters
%Tensor of parameters can be made by tt_addvecttobl in a simple way

%Write two new functions, done! :)


e=tt_tensor(e);
dm=tt_tensor(dm);
par=[];
for i=1:M
   par=kron(par,e);
end
for i=1:M
   w=[];
   for k=1:M
     if ( i == k )
        w=kron(w,dm);
     else
        w=kron(w,e);
     end
   end
   par=[par;w];
end
par=round(par,1e-9);


rhs=ones(size(mat{1},1),1);
rhs0=tt_ones(M*d0,2);
rhs0=tt_mkron({rhs,rhs0});
rhs0=tt_tensor(rhs0);


eps=1e-5;
%Simple iterative solver
tic;
x=solve_parametric3(mat,par,rhs0,[],1e-5,9);

toc;
return
x=tt_tensor(x);

%Compute the mean value in the centrum

ww2=tt_tensor(reshape(ww1,2*ones(1,d0)),1e-8);
ww3=[];
for i=1:M
  ww3=kron(ww3,ww2);
end
%Vast stupidity but fast programming
mean_m=zeros(n0);
for i=1:n0
for j=1:n0
    e0=zeros(n0); e0(i,j)=1; e0=e0(:); e0=tt_tensor(e0);
   z1=kron(e0,ww3);
   mean_m(i,j)=dot(x,z1);
    %e0=zeros(n0,n0); e0(n0/2,n0/2)=1; e0=e0(:); e0=tt_tensor(e0);
%z1=kron(e0,ww3);
%ww4=dot(x,z1);
end
end

fprintf('Mean value in the mid point is: %d \n',ww4);
%e0=tt_ones(x.d,x.n); e0=tt_tensor(e0);e0=e0/prod(x.n+1); %HEhe
%fprintf('NORM=%10.10f\n',dot(x,x)/(np+1)^M);
return
