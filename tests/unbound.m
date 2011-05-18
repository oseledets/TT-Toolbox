%Generate a simple algebraic unbounded grid 
% L=0.5;
% %The map is x/(x+L);
% d=10; n=2^d;
% h=1/n;
% %z=(0:n-1)*h;
% %x=L*z./(1-z);
% L=10;
% a=-20;
% b=20;
% h=(b-a)/(n-1);
% x=(0:n-1)*h;
% %max(x)
% [xx,yy]=meshgrid(x,x);
% ff=exp(-0.5*(xx-yy).^2);
% ff=exp(-0.5*xx.^2).*ff.*exp(-0.5*yy.^2);
% ff=reshape(ff,2*ones(1,2*d));
% p1=tt_matrix(ff,1e-8);
% p2=tt_tensor(ff,1e-8);
% p3=tt_tensor(reshape(ff,n,n),1e-8);

%Generate a simple "weighted" test
d=10; n=2^d;
a=-10;
b=10;
h=(b-a)/(n-1);
t=tt_x(d,2);t=tt_tensor(t);
e=tt_ones(d,2);e=tt_tensor(e);
t=a*e+t*h;
%x=kron(kron(t,e),e); y=kron(kron(e,t),e); z = kron(kron(e,e), t);
x=kron(t,e); y=kron(e,t);

p=(x+y).^2; p=round(p,1e-8);
%p=(x+y+z).^2;
%p1=x.^2+y.^2+z.^2; p1=round(p1,1e-8);

p1=x.^2+y.^2; p1=round(p1,1e-8);
%Generate the weight
w=funcross(p1,@(x) exp(-x),1e-8,p1,10);
p=funcross(p,@(x) exp(-x),1e-8,p,10);
p=x.^2+y.^2 + p; p =round(p,1e-8);
% keyboard;
p2=w.*p; 
p2=round(p2,1e-8);

% !!! Buggy place !!!! !!! Buggy place !!!! !!! Buggy place !!!!
% Use "full" carefully and only for small sizes!
% Otherwise - OOM killer will do you.
pf2=full(p2,[n^2, 1]);
w1=full(w,[n^2, 1]);
pf2=pf2./w1;

norm_pf2 = norm(pf2-full(p, [n^2,1]))/norm(full(p, [n^2,1]))



% close all;
% mesh(pf2-full(p,[n,n]));
% figure;
% mesh(full(round(p,1e-5),[n,n])-full(p,[n,n]));
%Generate the linear system and solve it
%rhs=diag(w)*p; %rhs=round(rhs,1e-8);
%sol=dmrg_solve2(diag(w),rhs,p,1e-5,1e-5,1000,10,[]);

