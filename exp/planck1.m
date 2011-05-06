d=8; n=2^d;
a=-12;b=12;
h=(b-a)/(n+1);
lp=tt_qlaplace_dd([d,d]);
lp=tt_matrix(lp); 
lp=lp/h^2;
z=tt_shf(d); z=tt_matrix(z);
e=tt_eye(2,d); e=tt_matrix(e);
c=(z-z')/(2*h);
cx=kron(c,e); cy=kron(e,c);
t=tt_x(d,2); t=tt_tensor(t);
e1=tt_ones(d,2); e1=tt_tensor(e1); 
t=t+e1; t=t*h; x0=a*e1+t; x0=round(x0,1e-10);
x=kron(x0,e1); y=kron(e1,x0);
dx=2*x; 
dy=2*(x.*y); %For the gaussian case
rs=x.^2+y.^2; rs=round(rs,1e-10);

%r1=funcrs2(rs,@(p) 1.0./sqrt(p),1e-13,rs,15); r1=round(r1,1e-10);
%dx=x.*r1; dx=round(dx,1e-8);
%dy=y.*r1; dy=round(dy,1e-8);

dx=diag(dx);
dy=diag(dy);



mat=-lp+cx*dx+cy*dy;



%norm(pp+pp')/norm(pp)
%return
mat=round(mat,1e-8);
%Explicit solution is known
st=tt_random(size(rs),ndims(rs),8); st=tt_tensor(st);
f=funcrs(rs,@(p) exp(-p),1e-13,st,10);
norm(mat*f)/norm(f)
%return
%f=funcrs2(rs,@(p) exp(-sqrt(p)),1e-13,rs,15);

%Create rhs;
%p0=cx*dx+cy*dy; p0=round(p0,1e-8);
%rhs=tt_mvk2(core(p0),core(f),1e-8);
%rhs=tt_tensor(rhs); rhs=-rhs;

%sol=tt_random(2,2*d,10); sol=tt_tensor(sol);
%rhs=tt_random(2,ndims(rhs),2);
%srhs=tt_tensor(rhs); rhs1=rhs; rhs=mat*rhs; rhs=round(rhs,1e-10);
tau=5e-2;
S=kron(e,e)-tau*mat; S=round(S,1e-8);
e2=kron(e1,e1);
rr=(x-5*e2).^2+(y-3*e2).^2;
rr=round(rr,1e-8);
sol=funcrs(rr,@(p) exp(-p),1e-13,rr,10); 
%rhs=tt_random(2,ndims(rhs
close all;
for i=1:60
    sol1=reshape(full(sol),n,n); 
%mesh(sol1);
contourf(sol1);
drawnow;
    rhs0=sol;
%    sol=dmrg_solve2(core(S),core(rhs0),core(sol),1e-7,1e-7,100,8);
sol=dmrg_solve2(S,rhs0,sol,1e-7,1e-7,100,8);
sol=tt_tensor(sol);
%norm(S*sol-rhs0)/norm(rhs0)
%norm(sol)
%sol1=reshape(full(sol),n,n); 
%mesh(sol1);
%contourf(sol1);
%drawnow;
%norm(mat*sol)/norm(sol)
%dot(sol,f)/(norm(sol)*norm(f))
%keyboard;
%keyboard;
end
