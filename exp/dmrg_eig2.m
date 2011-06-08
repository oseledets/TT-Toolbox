function [x]=dmrg_eig2(a,x0,eps,nswp)
%[X]=DMRG_EIG2(A,X0,EPS)
%Minimal eigenvalue of a QTT matrix
%It uses "Density Matrix Normalization Group" Ideology in the following
%way
%It mimics Krylov matrix-by-vector product
if (nargin<4)||(isempty(nswp))
    nswp=50;
end;
rmin=1; 
kickrank=10;
verb=false;
%eps=1e-4;
rmax=100;
%verb=true;
d=size(a,1);
x=cell(d,1);
sz=zeros(d,1);
%Init solution 
% for i=1:d
%   sz(i) = size(a{i},1);
%   x{i}=randn(sz(i),1);
% end
%x=tt_ones(d,2);
x=x0;
%x=tt_scal(y,5.0); %Initialize to right-hand side
%x0=x;
eps0=eps/sqrt(d-1);
phx=cell(d,1); %Phi matrices for X^T A X
swp=1;
not_converged = true;
while ( swp < nswp && not_converged )
%      not_converged = false;
     er_max=0;
 %x0=tt_random(tt_size(x),size(x,1),1);
 %   x=tt_add(x,x0);
  
    %   z=tt_mv(a,x);
 %   er=tt_dist2(z,y);
 %   nrm=sqrt(tt_dot(z,z));
 %   fprintf('er=%3.2e\tnrm=%3.1e\trel_er=%3.2e\n',er,nrm,er/nrm);
%Sweep starts from right-to-left orthogonalization of x
%Compute right-to-left qr and phi matrices
mat=x{d};
[q,rv]=qr(mat,0);
x{d}=q;
%A(n,m,ra))*x(m,rxm)*x(n,rxn1) -> phx(ra,rxm,rxn)
%y(m,ry)*x(m,rx) -> phy(ry,rx)
core=a{d}; n=size(core,1); m=size(core,2); ra=size(core,3);
mat_x=x{d}; rx=size(mat_x,2);
core=permute(core,[1,3,2]); 
core=reshape(core,[n*ra,m]);
core=core*mat_x; %is n*ra*rx1
core=reshape(core,[n,ra,rx]); %n*ra*rxm
core=permute(core,[2,3,1]); %ra*rxm*n
core=reshape(core,[ra*rx,n]);
phx{d}=core*mat_x; % ra*rxm*rxn
phx{d}=reshape(phx{d},[ra,rx,rx]);

%Back reorthogonalization
for i=(d-1):-1:2
  
    core=x{i};
    core=ten_conv(core,3,rv');
    ncur=size(core,1);
    r2=size(core,2);
    r3=size(core,3);
    core=permute(core,[1,3,2]);
    core=reshape(core,[ncur*r3,r2]);
    [x{i},rv]=qr(core,0);
    rnew=min(r2,ncur*r3);
    x{i}=reshape(x{i},[ncur,r3,rnew]);
    x{i}=permute(x{i},[1,3,2]);
    %A(n,m,ra1,ra2)*x(m,rxm1,rxm2)*x(n,rxn1,rxn2)*phx(ra2,rxm2,rxn2)->
    %phx(ra1,rxm1,rxn1)
    %Decipher sizes
    core=a{i}; n=size(core,1); m=size(core,2); ra1=size(core,3); ra2=size(core,4);
    core_x=x{i}; 
    rx3=size(core_x,2); rx1=rx3; rx4=size(core_x,3); rx2=rx4;
    ph=phx{i+1};
    core_x=reshape(core_x,[n*rx3,rx4]); 
    ph=reshape(ph,[ra2*rx2,rx4]);
    ph=core_x*ph'; %ph is n*rx3*ra2*rx2
    ph=reshape(ph,[n,rx3,ra2,rx2]);
    %A(n,m,ra1,ra2)*ph(n,rx3,ra2,rx2) over (n,ra2)
    core=permute(core,[2,3,1,4]); core=reshape(core,[m*ra1,n*ra2]);
    ph=permute(ph,[1,3,2,4]); ph=reshape(ph,[n*ra2,rx3*rx2]);
    ph=core*ph; %ph is m*ra1*rx3*rx2, convolve with x(m,rx1,rx2)
    ph=reshape(ph,[m,ra1,rx3,rx2]); 
    ph=permute(ph,[2,3,1,4]);
    ph=reshape(ph,[ra1*rx3,m*rx2]);
    core_x=permute(x{i},[2,1,3]); 
    core_x=reshape(core_x,[rx1,m*rx2]);
    ph=ph*core_x'; %ph is ra1*rx3*rx1
    ph=reshape(ph,[ra1,rx3,rx1]);
    ph=permute(ph,[1,3,2]); 
    phx{i}=ph;
    %x(n,rx3,rx4)*y(n,ry3,ry4)*phy(ry4,rx4)
   end
%Now start main dmrg iteration
%First, the first core
%A1(n1,m1,ra1)*A2(n2,m2,ra1,ra2)*phx(ra2,rxm2,rxn2)  --- 
%matrix of size B(rxn2*n1n2,m1m2*rxm2)
%Unknown: X1(m1,rx1)*X2(m2,rx1,rx2) -> m1*m2*rxm2
%Decipher sizes
a1=a{1}; n1=size(a1,1); m1=size(a1,2); ra1=size(a1,3);
a2=a{2}; n2=size(a2,1); m2=size(a2,2); ra2=size(a2,4);
phx2=phx{3};
rxm2=size(phx2,2); 
rxn2=size(phx2,3);
a1=reshape(a1,[n1*m1,ra1]);
a2=permute(a2,[3,1,2,4]);
a2=reshape(a2,[ra1,n2*m2*ra2]);
b=a1*a2;
b=reshape(b,[n1*m1*n2*m2,ra2]);
phx2=reshape(phx2,[ra2,rxm2*rxn2]);
b=b*phx2;
b=reshape(b,[n1,m1,n2,m2,rxm2,rxn2]);
b=permute(b,[1,3,6,2,4,5]);
b=reshape(b,[n1*n2*rxn2,m1*m2*rxm2]);

% 
%  %Here we also have a previous solution
%   x1=x{1};
%   x2=x{2}; 
%   rxm1=size(x1,2);
%   x2=permute(x2,[2,1,3]); x2=reshape(x2,[rxm1,n2*rxm2]);
%   sol_prev=x1*x2; 
%   sol_prev=reshape(sol_prev,[n1*n2*rxm2,1]);
%   
%  
   
   [sol,ev]=eigs((b+b')*0.5,1,'sm');

   
%   
%   
%   er=acos(abs(dot(sol_prev,sol))/(norm(sol)*norm(sol_prev)));
%   er=abs(er);
%   if ( er > eps )
%     not_converged=true;
%     er_max = max(er_max,er);
%   end
   
%keyboard;

%sol=gmres(b,rhs,10,1e-6,1000);



%tt_dist2(tt_mv(a,x),y)
%sol is m1*m2*rxm2
sol=reshape(sol,[m1,m2*rxm2]);
[u,s,v]=svd(sol,'econ');
s=diag(s);
nrm=norm(s);
rold=size(x{1},2);
r=my_chop2(s,eps0*nrm);
r=max([rmin,r]);
r=min(r,rmax);
r=min(r,size(s,1));

if ( verb ) 
 fprintf('We can push rank %d to %d \n',1,r);
end
u=u(:,1:r);
x{1}=u;
 v=v(:,1:r)*diag(s(1:r));
 v=reshape(v,[m2,rxm2,r]);
 v=permute(v,[1,3,2]);
 x{2}=v;
 %keyboard;
%Recalculate phx{1},phy{1};
%A(n,m,ra)*X(m,rxm)*X(n,rxn)-> phx(rxm,rxn,ra)
%Y(n,ry)*X(n,rxn) -> phy(rxn,ry)
%Decipher sizes
core=a{1}; n=size(core,1); m=size(core,2); ra=size(core,3);
matx=x{1}; rxm=size(matx,2); rxn=size(matx,2);
core=reshape(core,[n,m*ra]);
ph=matx'*core; %ph is rxn*m*ra
ph=reshape(ph,[rxn,m,ra]);
ph=permute(ph,[2,1,3]);
ph=reshape(ph,[m,rxn*ra]);
ph=matx'*ph;
phx{1}=reshape(ph,[rxm,rxn,ra]);

for i=2:d-2
%% For main iteration:
 
%ph1(rxm1,rxn1,ra1)*A1(n1,m1,ra1,ra2)*A2(n2,m2,ra2,ra3)*ph2(ra3,rxm3,rxn3)
 %is a matrix of size (rxn1*rxn2) n1n2 m1m2 (rxm1*rxm3)
 %right-hand side:
 %phy1(rxn1,ry1)*Y1(n1,ry1,ry2)*Y2(n2,ry2,ry3)*phy2(ry3,rxn3)
 %has size rxn1*n1n2*rxn3
 
 %Decipher sizes
 a1=a{i}; n1=size(a1,1); m1=size(a1,2); ra1=size(a1,3); ra2=size(a1,4);
 a2=a{i+1}; n2=size(a2,1); m2=size(a2,2); ra3=size(a2,4);
 ph1=phx{i-1}; rxm1=size(ph1,1); rxn1=size(ph1,2);
 ph3=phx{i+2}; rxm3=size(ph3,2); rxn3=size(ph3,3);
 
 
 %We had a previous solution:
 %X1(m1,rx1,rx2)*X2(m2,rx2,rx3)
 x1=x{i}; x2=x{i+1}; 
 rx1=size(x1,2); rx2=size(x1,3); rx3=size(x2,3);
 x2=permute(x2,[2,1,3]); x2=reshape(x2,[rx2,m2*rx3]);
 x1=permute(x1,[2,1,3]);
 x1=reshape(x1,[rx1*m1,rx2]);
 sol_prev=x1*x2; 
 sol_prev=reshape(sol_prev,[rxm1*m1*m2*rxm3,1]);
 
 
 
%ph3 is ra2*n2*m2*rxm3*rxn3
%  ph1=reshape(ph1,[rxm1*rxn1*n1*m1,ra2]);
%  ph3=reshape(ph3,[ra2,n2*m2*rxm3*rxn3]);
%  b=ph1*ph3; 
%  %b is rxm1*rxn1*n1*m1*n2*m2*rxm3*rxn3
%  b=reshape(b,[rxm1,rxn1,n1,m1,n2,m2,rxm3,rxn3]);
%  % is (rxn1*rxn3)*n1n2m1m2 (rxm1 rxm3)
%  b=permute(b,[2,3,5,8,1,4,6,7]);
%  b=reshape(b,[rxn1*n1*n2*rxn3,rxm1*m1*m2*rxm3]);
 
 % opts.issym = 1;
 % opst.v0=sol_prev;
 %[sol,ev]=eigs(@(x)bfun1(x,a1,a2,ph1,ph3,n1,m1,n2,m2,ra1,ra2,ra3,rxn1,rxm1,rxn3,rxm3),rxn1*n1*n2*rxn3,1,'SM',opts);
 if ( rxn1*n1*n2*rxn3 <= 6)
 
 %Decipher sizes
 a1=a{i}; n1=size(a1,1); m1=size(a1,2); ra1=size(a1,3); ra2=size(a1,4);
 a2=a{i+1}; n2=size(a2,1); m2=size(a2,2); ra3=size(a2,4);
 ph1=phx{i-1}; rxm1=size(ph1,1); rxn1=size(ph1,2);
 ph3=phx{i+2}; rxm3=size(ph3,2); rxn3=size(ph3,3);
 a1=permute(a1,[3,1,2,4]); a1=reshape(a1,[ra1,n1*m1*ra2]);
 ph1=reshape(ph1,[rxm1*rxn1,ra1]);
 ph1=ph1*a1; %ph1 is rxm1*rxn1*n1*m1*ra2;
 ph_save=ph1;
 a2=permute(a2,[3,1,2,4]); a2=reshape(a2,[ra2*n2*m2,ra3]); 
 ph3=reshape(ph3,[ra3,rxm3*rxn3]);
 ph3=a2*ph3; %ph3 is ra2*n2*m2*rxm3*rxn3
 ph1=reshape(ph1,[rxm1*rxn1*n1*m1,ra2]);
 ph3=reshape(ph3,[ra2,n2*m2*rxm3*rxn3]);
 b=ph1*ph3; 
 %b is rxm1*rxn1*n1*m1*n2*m2*rxm3*rxn3
 b=reshape(b,[rxm1,rxn1,n1,m1,n2,m2,rxm3,rxn3]);
 % is (rxn1*rxn3)*n1n2m1m2 (rxm1 rxm3)
 b=permute(b,[2,3,5,8,1,4,6,7]);
 b=reshape(b,[rxn1*n1*n2*rxn3,rxm1*m1*m2*rxm3]);
 [sol,ev]=eigs((b+b')*0.5,1,'sm');
 else
 [sol,ev]=lobpcg(sol_prev,@(x)bfun1(x,a1,a2,ph1,ph3,n1,m1,n2,m2,ra1,ra2,ra3,rxn1,rxm1,rxn3,rxm3));
  fprintf('ev=%4.10f \n',ev);
 
 end
 %keyboard;
 %fprintf('ev=%3.4f \n',ev);
 
 
 % [sol,ev]=eigs(b,1,'sm');
  %er=norm(sol_prev-sol)/norm(sol);
  er=acos(abs(dot(sol_prev,sol))/(norm(sol)*norm(sol_prev)));
  er=abs(er);
  if ( er > eps )
    not_converged=true;
    er_max = max(er_max,er);
  end
 
%   fprintf('block %d, norm_imag=%3.3e, b-bT=%3.3e, er=%3.3e\n', i, norm(imag(sol)), norm(b-b'), er);

 
 
 
 %sol is rxm1*m1*m2*rxm3
 sol=reshape(sol,[rxm1*m1,m2*rxm3]);
 [u,s,v]=svd(sol,'econ');
  s=diag(s);  
  %fprintf('norm=%18f \n',norm(s));
  %fprintf('tensor norm=%3.2e \n',norm(s));
  rold=size(x{i},3);
  r=my_chop2(s,eps0*norm(s));
  r=max([r,rmin]);
  r=min(r,rmax);
  r=min(r,size(s,1));
  %if ( 
  if ( verb )
  fprintf('We can push rank %d to %d \n',i,r);
  end
  
  % random kick
  u = u(:,1:r);
  u=[u, randn(size(u,1), kickrank)*max(abs(u(:)))];
  v=v(:,1:r)*diag(s(1:r));
  % random kick
  v = [v, zeros(size(v,1), kickrank)];
  [u,rv]=qr(u,0);
  r = size(u,2);
  v = v*(rv.');
  
  u=reshape(u,[rxm1,m1,r]);
  x{i}=permute(u,[2,1,3]);
  %Recompute phx and phy
  %phx(rxm2,rxn2,ra2) =
  %ph1(rxm1,rxn1,ra1)*A1(n1,m1,ra1,ra2)
  %X(m1,rxm1,rxm2)*X(n1,rxn1,rxn2) -> ph(rxm2,rxn2,ra2)
  
  %Recalculate everything from scratch to save memory
  ph1=phx{i-1};
  
  x1=x{i};
  rxm2=size(x1,3);
  x1=permute(x1,[2,1,3]);
  ph1=reshape(ph1,[rxm1,rxn1*ra1]);
  x1=reshape(x1,[rxm1,m1*rxm2]);
  ph1=ph1'*x1; %ph1 is now rxn1*ra1*m1*rxm2
  %Sum with m1,ra1
  ph1=reshape(ph1,[rxn1,ra1,m1,rxm2]);
  ph1=permute(ph1,[2,3,1,4]);
  ph1=reshape(ph1,[ra1*m1,rxn1*rxm2]);
  a1=a{i};
  a1=permute(a1,[3,2,1,4]);
  a1=reshape(a1,[ra1*m1,n1*ra2]);
  ph1=ph1'*a1; %is rxn1*rxm2*n1*ra2
  x1=x{i}; rxn2=size(x1,3); x1=reshape(x1,[n1*rxn1,rxn2]);
  %Over n1,rxn1
  ph1=reshape(ph1,[rxn1,rxm2,n1,ra2]);
  ph1=permute(ph1,[2,4,3,1]);
  ph1=reshape(ph1,[rxm2*ra2,n1*rxn1]);
  ph1=ph1*x1; %is rxm2*ra2*rxn2;
  ph1=reshape(ph1,[rxm2,ra2,rxn2]);
  ph1=permute(ph1,[1,3,2]);
  phx{i}=ph1;
  %phy(rxn1,ry1)*X(n1,rxn1,rxn2)*Y(n1,ry1,ry2)-> phy(rxn2,ry2)
  v=reshape(v,[m2,rxm3,r]);
  v=permute(v,[1,3,2]);
  x{i+1}=v;
 % tt_dist2(tt_mv(a,x),y)
 % tt_dist2(tt_mv(a,x0),y)
 % keyboard;
end
  

  %And compute the last core 
  %ph1(rxm1,rxn1,ra1)*a1(n1,m1,ra1,ra2)*a2(n2,m2,ra2)
  %Decipher sizes
  a1=a{d-1}; n1=size(a1,1); m1=size(a1,2); ra1=size(a1,3); ra2=size(a1,4);
  a2=a{d}; n2=size(a2,1); m2=size(a2,2);
  ph1=phx{d-2}; 
  rxm1=size(ph1,1); rxn1=size(ph1,2);
  a1=permute(a1,[3,1,2,4]); a1=reshape(a1,[ra1,n1*m1*ra2]);
  ph1=reshape(ph1,[rxm1*rxn1,ra1]);
  ph1=ph1*a1; %ph1 is rxm1*rxn1*n1*m1*ra2
  ph1=reshape(ph1,[rxm1*rxn1*n1*m1,ra2]);
  a2=permute(a2,[3,1,2]); a2=reshape(a2,[ra2,n2*m2]);
  ph1=ph1*a2; 
  %ph1 is rxm1*rxn1*n1*m1*n2*m2
  ph1=reshape(ph1,[rxm1,rxn1,n1,m1,n2,m2]);
  ph1=permute(ph1,[3,5,2,4,6,1]);
  b=reshape(ph1,[n1*n2*rxn1,m1*m2*rxm1]);
  
  %Here we also have a previous solution
  x1=x{d-1};
  x2=x{d}; 
  rxn2=size(x2,2);
  sol_prev=reshape(x1,[n1*rxn1,rxn2])*x2'; 
  sol_prev=reshape(sol_prev,[n1,rxn1,n2]);
  sol_prev=permute(sol_prev,[1,3,2]);
  sol_prev=reshape(sol_prev,[n1*n2*rxn1,1]);
  
 
  
  [sol,ev]=eigs((b+b')*0.5,1,'sm');
  
  
  
  er=acos(abs(dot(sol_prev,sol))/(norm(sol)*norm(sol_prev)));
  er=abs(er);
  if ( er > eps )
    not_converged=true;
    er_max = max(er_max,er);
  end
   
  
  %sol is m1*m2*rxm1
  sol=reshape(sol,[m1,m2,rxm1]);
  sol=permute(sol,[1,3,2]);
  sol=reshape(sol,[m1*rxm1,m2]);
  [u,s,v]=svd(sol,'econ');
  s=diag(s);  
  rold=size(x{d-1},3);
  r=my_chop2(s,eps0*norm(s));
  r=max([r,rmin]);
  r=min(r,rmax);
  r=min(r,size(s,1));
  %if ( 
  if ( verb )
  fprintf('We can push rank %d to %d \n',1,r);
  end
  u=u(:,1:r);
  u=reshape(u,[m1,rxm1,r]);
  x{d-1}=u;
  v=v(:,1:r)*diag(s(1:r));
  x{d}=v;
  fprintf('swp=%d er_max=%3.2e \n',swp,er_max);
  swp=swp+1;
  
end
%fprintf('Last stand \n');
%keyboard;
return
end

function [y]=bfun1(x,a1,a2,ph1,ph2,n1,m1,n2,m2,ra1,ra2,ra3,rxn1,rxm1,rxn3,rxm3)
%[Y]=BFUN1(x,a1,a2,ph1,ph2,n1,m1,n2,m2,ra1,ra2,ra3,rxn1,rxm1,rxn2,rxm2)
%This MEGA function computes a single matrix-by-vector product
%ph1(rxm1,rxn1,ra1)*A1(n1,m1,ra1,ra2)*A2(n2,m2,ra2,ra3)*ph2(ra3,rxm3,rxn3)
%X(rxm1,m1,m2,rxm3)

y=reshape(ph1,[rxm1,rxn1*ra1]);
x=reshape(x,[rxm1,m1*m2*rxm3]);
y=y'*x; %Is rxn1*ra1*m1*m2*rxm3
%Convolve over m1,ra1
a1=permute(a1,[2,3,1,4]); %is [m1*ra1,n1*ra2]);
y=reshape(y,[rxn1,ra1,m1,m2,rxm3]);
y=permute(y,[3,2,1,4,5]);
y=reshape(y,[m1*ra1,rxn1*m2*rxm3]);
a1=reshape(a1,[m1*ra1,n1*ra2]);
y=y'*a1; %is rxn1*m2*rxm3*n1*ra2, over m2*ra2 with a2
y=reshape(y,[rxn1,m2,rxm3,n1,ra2]);
y=permute(y,[1,4,3,2,5]);
y=reshape(y,[rxn1*n1*rxm3,m2*ra2]);
a2=permute(a2,[2,3,1,4]);
a2=reshape(a2,[m2*ra2,n2*ra3]);
y=y*a2; %y is rxn1*n1*rxm3*n2*ra3
y=reshape(y,[rxn1*n1,rxm3,n2,ra3]);
y=permute(y,[1,3,4,2]);
y=reshape(y,[rxn1*n1*n2,ra3*rxm3]);
ph2=reshape(ph2,[ra3*rxm3,rxn3]);
y=y*ph2;
y=reshape(y,[rxn1*n1*n2*rxn3,1]);
return
end