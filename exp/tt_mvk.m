function [y]=tt_mvk(a,x,eps)
%[A]=TT_MVK(A,X,EPS)
%Matrix-by-vector product by Krylov-type Block algorithm,
%which tries to compute everything fast, by increasing rank as far as
%possible
rmin=1; nswp=6;
radd=0;
verb=false;
%verb=true;
d=size(a,1);
y=cell(d,1);
sz=zeros(d,1);
%Init solution 
for i=1:d
  sz(i) = size(a{i},1);
  y{i}=randn(sz(i),1);
end
%First initialization is done as one sweep of tt_mvrank, left-right and
%right-left computes "correct" y (with normalization, also)
[y,ph]=tt_mvr(a,x,y);
%Now start main iteration (computing b-matrices and their orthogonalization
%etc)

%ph will be rx * ra * ry 
for swp=1:nswp

%First, the first mode
%A(n1,m1,ra1)*X(m1,rx1)*A(n2,m2,ra1,ra2)*X(m2,rx1,rx2)*ph(rx2,ra2,ry2)
%Output will be B(n1,n2,ry2), will separate n1; n2ry2
%Decipher sizes
core1=a{1}; n1=size(core1,1); m1=size(core1,2); ra1=size(core1,3);
core2=a{2}; n2=size(core2,2); m2=size(core2,2); ra2=size(core2,4);
x1=x{1}; rx1=size(x1,2); 
x2=x{2}; rx2=size(x2,3);
phi2=ph{3};
ry2=size(phi2,3);
%X(m2,rx1,rx2)*ph(rx2,ra2,ry2);
phi2=reshape(x2,[m2*rx1,rx2])*reshape(phi2,[rx2,ra2*ry2]);
phi2=reshape(phi2,[m2,rx1,ra2,ry2]);
phi2=permute(phi2,[1,3,2,4]);
core2=permute(core2,[1,3,2,4]);
phi2=reshape(core2,[n2*ra1,m2*ra2])*reshape(phi2,[m2*ra2,rx1*ry2]);
phi2=reshape(phi2,[n2,ra1,rx1,ry2]);
phi2=permute(phi2,[2,3,1,4]);
phi2=reshape(phi2,[ra1*rx1,n2*ry2]);
core1=permute(core1,[1,3,2]);
x1=reshape(core1,[n1*ra1,m1])*x1;
mat=reshape(x1,[n1,ra1*rx1])*phi2;
mat=reshape(mat,[n1,n2*ry2]);
%Project mat into the right bases

%mat is n1*n2*ry2
%y0=y{1}; r0=size(y0,2);
%y1=y{3}; 
%ry3=size(y1,3);
%y1 is [n2,ry2,ry3] and  is krivo-orthogonal
%y1=permute(y1,[1,3,2]); y1=reshape(y1,[n2*ry3,ry2]);
%e0=eye(n1); e1=eye(n2*ry2);
%keyboard;
nrm=norm(mat,'fro'); %Estimate of the tensor norm
%phi=y0'*mat*y1;
%mat1=mat-y0*phi*y1';
%mat=(e0-y0*y0')*mat*(e1-y1*y1');
%Determine eps-rank of mat
%keyboard;
[u,s,v]=svd(mat);
s=diag(s);
rold=size(y{1},2);
r=my_chop2(s,eps*nrm);
r=max([rmin,r,rold+radd]);
r=min(r,size(s,1));
if ( verb ) 
 fprintf('We can push rank %d to %d \n',1,r);
end
u=u(:,1:r);
v=v(:,1:r);
v=v*diag(s(1:r));
y{1}=u;
%y{2}=v; %v is n2*ry2*ry1
%v=reshape(v,[n2,ry2,r]);
%v=permute(v,[1,3,2]); v=reshape(v,[n2*r,ry2]); [v,dmp]=qr(v,0);
%v=reshape(v,[n2,r,ry2]);
%y{2}=v;
%Recalculate ph{1};
%ph{1}=A(n,m,ra)*X(m,rx)*Y(n,ry) -> (ry,ra,rx)
%Decipher sizes
core=a{1}; n=size(core,1); m=size(core,2); ra=size(core,3);
x1=x{1}; rx=size(x1,2);
y1=y{1}; ry=size(y1,2);
core=permute(core,[1,3,2]); core=reshape(core,[n*ra,m]);
phi1=core*x1; %phi is n*ra*rx
phi1=y1'*reshape(phi1,[n,ra*rx]);
ph{1}=reshape(phi1,[ry,ra,rx]);
 for i=2:d-2
     
  %That will be not very easy
  %First matrix B is computed
  %ph1(ry1,ra1,rx1)*A1(n1,m1,ra1,ra2)*X1(m1,rx1,rx2)
  %A2(n2,m2,ra2,ra3)*X2(m2,rx2,rx3)*ph2(rx3,ra3,ry3)
  %Decipher sizes
  core1=a{i};  n1=size(core1,1); m1=size(core1,2); ra1=size(core1,3); ra2=size(core1,4);
  core2=a{i+1}; n2=size(core2,1); m2=size(core2,2); ra3=size(core2,4);
  x1=x{i}; rx1=size(x1,2); rx2=size(x1,3);
  x2=x{i+1}; rx3=size(x2,3); 
  ph1=ph{i-1}; ry1=size(ph1,1);
  ph2=ph{i+2}; ry3=size(ph2,3);
  %ph1(ry1,ra1,rx1)*X1(m1,rx1,rx2)
  x1=permute(x1,[2,1,3]); x1=reshape(x1,[rx1,m1*rx2]);
  ph1=reshape(ph1,[ry1*ra1,rx1]);
  ph1=ph1*x1; %ph1 is ry1*ra1*m1*rx2
  %ph1(ry1,ra1,m1,rx2)*A1(n1,m1,ra1,ra2) m1,ra1
  ph1=reshape(ph1,[ry1,ra1,m1,rx2]);
  ph1=permute(ph1,[4,1,3,2]);
  ph1=reshape(ph1,[rx2*ry1,m1*ra1]);
  core1=permute(core1,[2,3,4,1]);
  core1=reshape(core1,[m1*ra1,ra2*n1]);
  ph1=ph1*core1; 
  ph1=reshape(ph1,[rx2,ry1,ra2,n1]);
  %ph1 is rx2*ry1*ra2*n1 %This one will be needed afterwards
  ph_save=ph1;
  %The same goes for second
  ph2=reshape(x2,[m2*rx2,rx3])*reshape(ph2,[rx3,ra3*ry3]);
  %ph2 is m2*rx2*ra3*ry3
  ph2=reshape(ph2,[m2,rx2,ra3,ry3]);
  %A2(n2,m2,ra2,ra3)*ph2(m2,rx2,ra3,ry3) over m2,ra3
  core2=permute(core2,[1,3,2,4]); core2=reshape(core2,[n2*ra2,m2*ra3]);
  ph2=permute(ph2,[1,3,2,4]); ph2=reshape(ph2,[m2*ra3,rx2*ry3]);
  ph2=core2*ph2;
  %ph2 is n2*ra2*rx2*ry3
  ph2=reshape(ph2,[n2,ra2,rx2,ry3]);
  %Recall, ph1 is rx2*ry1*ra2*n1 convolution is overs rx2,ra2
  %Matrix B is ry1*n1*n2*ry2
  ph1=permute(ph1,[2,4,3,1]); %Now ph1 is (ry1*n1)*(ra2*rx2)
  ph1=reshape(ph1,[ry1*n1,ra2*rx2]);
  ph2=permute(ph2,[1,4,2,3]);
  ph2=reshape(ph2,[n2*ry3,ra2*rx2]);
  mat=ph1*ph2';
  [u,s,v]=svd(mat,'econ');
   s=diag(s);
  
%fprintf('norm=%18f \n',norm(s));
  %fprintf('tensor norm=%3.2e \n',norm(s));
  rold=size(y{i},3);
  r=my_chop2(s,eps*norm(s));
r=max([r,rmin,rold+radd]);
r=min(r,size(s,1));
  %if ( 
  if ( verb )
  fprintf('We can push rank %d to %d \n',i,r);
  end
  u=u(:,1:r);
  %v=v(:,1:r)*diag(s(1:r));
  u=reshape(u,[ry1,n1,r]);
  y{i}=permute(u,[2,1,3]);
  %v=reshape(v,[n2,ry3,r]);
  %y{i+1}=permute(v,[1,3,2]);
  %Finally, need to recompute ph{i}
  %will have size ry
  %A(n1,m1,ra1,ra2)*X(m1,rx1,rx2)*Y(n1,ry1,ry2)*ph{i-1}(ry1,ra1,rx1)
  %ph_save(rx2,ry1,ra2,n1)*Y(n1,ry1,ry2) 
  ph_save=reshape(ph_save,[rx2,ry1,ra2,n1]);
  ph_save=permute(ph_save,[4,2,3,1]);
  ph_save=reshape(ph_save,[n1*ry1,ra2*rx2]);  %HERE WE CAN HAVE A FAULT
  ph_save=reshape(y{i},[n1*ry1,r])'*ph_save;
  ph{i}=reshape(ph_save,[r,ra2,rx2]);
 end
 %Now the last mode 
 %ph1(ry1,ra1,rx1)*A1(n1,m1,ra1,ra2)*X1(m1,rx1,rx2)
 %                *A2(n2,m2,ra2)*X2(m2,rx2)
 %Decipher sizes
 core1=a{d-1}; n1=size(core1,1); m1=size(core1,2); ra1=size(core1,3); ra2=size(core1,4);
 core2=a{d}; n2=size(core2,1); m2=size(core2,2);
 x1=x{d-1}; rx1=size(x1,2); rx2=size(x1,3);
 x2=x{d};
 ph1=ph{d-2};
 ry1=size(ph1,1);
 %X1(m1,rx1,rx2)*ph1(ry1,ra1,rx1)
 %ph1(ry1,ra1,rx1)*X1(m1,rx1,rx2)
 x1=permute(x1,[2,1,3]); x1=reshape(x1,[rx1,m1*rx2]);
 ph1=reshape(ph1,[ry1*ra1,rx1]);
 ph1=ph1*x1; %ph1 is ry1*ra1*m1*rx2
 %ph1(ry1,ra1,m1,rx2)*A1(n1,m1,ra1,ra2) m1,ra1
 ph1=reshape(ph1,[ry1,ra1,m1,rx2]);
 ph1=permute(ph1,[4,1,3,2]);
 ph1=reshape(ph1,[rx2*ry1,m1*ra1]);
 core1=permute(core1,[2,3,4,1]);
 core1=reshape(core1,[m1*ra1,ra2*n1]);
 ph1=ph1*core1; 
 ph1=reshape(ph1,[rx2,ry1,ra2,n1]);
 %ph1 is rx2*ry1*ra2*n1 %This one will be needed afterwards
 %ph1=reshape
 %ph_save=ph1;
 
 %A2(n2,m2,ra2)*X2(m2,rx2) -> ph2(n2,ra2,x2)
 core2=permute(core2,[1,3,2]); core2=reshape(core2,[n2*ra2,m2]);
 ph2=core2*x2; %ph2 is n2*ra2*rx2
 ph_save=ph2; 
 ph2=reshape(ph2,[n2,ra2*rx2]);
 ph1=reshape(ph1,[rx2,ry1,ra2,n1]); %The needed result is ry1*n1*n2
 ph1=permute(ph1,[2,4,3,1]); ph1=reshape(ph1,[ry1*n1,ra2*rx2]);
 mat=ph1*ph2'; 
 [u,s,v]=svd(mat,'econ');
 s=diag(s);
 r=my_chop2(s,eps*norm(s));
 rold=size(y{d},2);
r=max([r,rmin,rold]);
r=min(r,size(s,1));
 %fprintf('We can push rank %d to %d \n',d-1,r);
 u=u(:,1:r);
 v=v(:,1:r); 
 %v=v;
 u=u*diag(s(1:r));
 u=reshape(u,[ry1,n1,r]);
 u=permute(u,[2,1,3]);
 y{d-1}=u;
 y{d}=v;
 %Recalculate ph{d}
 %A2(n2,m2,ra2)*X(m2,rx2)*Y(n2,ry2} %ph_save is n2*ra2*x2, result
 %rx2*ra2*ry2
 ph_save=reshape(ph_save,[n2,ra2,rx2]);
 ph_save=permute(ph_save,[1,3,2]); ph_save=reshape(ph_save,[n2,rx2*ra2]);
 ph_save=ph_save'*v;
 ph_save=reshape(ph_save,[rx2,ra2,r]);
 ph{d}=ph_save;
 %And start backwards
 for i=d-2:-1:2
     %Decipher sizes
  core1=a{i};  n1=size(core1,1); m1=size(core1,2); ra1=size(core1,3); ra2=size(core1,4);
  core2=a{i+1}; n2=size(core2,1); m2=size(core2,2); ra3=size(core2,4);
  x1=x{i}; rx1=size(x1,2); rx2=size(x1,3);
  x2=x{i+1}; rx3=size(x2,3); 
  ph1=ph{i-1}; ry1=size(ph1,1);
  ph2=ph{i+2}; ry3=size(ph2,3);
  %ph1(ry1,ra1,rx1)*X1(m1,rx1,rx2)
  x1=permute(x1,[2,1,3]); x1=reshape(x1,[rx1,m1*rx2]);
  ph1=reshape(ph1,[ry1*ra1,rx1]);
  ph1=ph1*x1; %ph1 is ry1*ra1*m1*rx2
  %ph1(ry1,ra1,m1,rx2)*A1(n1,m1,ra1,ra2) m1,ra1
  ph1=reshape(ph1,[ry1,ra1,m1,rx2]);
  ph1=permute(ph1,[4,1,3,2]);
  ph1=reshape(ph1,[rx2*ry1,m1*ra1]);
  core1=permute(core1,[2,3,4,1]);
  core1=reshape(core1,[m1*ra1,ra2*n1]);
  ph1=ph1*core1; 
  ph1=reshape(ph1,[rx2,ry1,ra2,n1]);
  %ph1 is rx2*ry1*ra2*n1 %This one will be needed afterwards
  %The same goes for second
  ph2=reshape(x2,[m2*rx2,rx3])*reshape(ph2,[rx3,ra3*ry3]);
  %ph2 is m2*rx2*ra3*ry3
  ph2=reshape(ph2,[m2,rx2,ra3,ry3]);
  %A2(n2,m2,ra2,ra3)*ph2(m2,rx2,ra3,ry3) over m2,ra3
  core2=permute(core2,[1,3,2,4]); core2=reshape(core2,[n2*ra2,m2*ra3]);
  ph2=permute(ph2,[1,3,2,4]); ph2=reshape(ph2,[m2*ra3,rx2*ry3]);
  ph2=core2*ph2;
  %ph2 is n2*ra2*rx2*ry3
  ph_save=ph2;
 
  ph2=reshape(ph2,[n2,ra2,rx2,ry3]);
  %Recall, ph1 is rx2*ry1*ra2*n1 convolution is overs rx2,ra2
  %Matrix B is ry1*n1*n2*ry3
  ph1=permute(ph1,[2,4,3,1]); %Now ph1 is (ry1*n1)*(ra2*rx2)
  ph1=reshape(ph1,[ry1*n1,ra2*rx2]);
  ph2=permute(ph2,[1,4,2,3]);
  ph2=reshape(ph2,[n2*ry3,ra2*rx2]);
  mat=ph1*ph2';
  [u,s,v]=svd(mat,'econ');
  s=diag(s);
  rold=size(y{i},3);
  r=my_chop2(s,eps*norm(s));
 r=max([r,rmin,rold+radd]); 
r=min(r,size(s,1));
  %fprintf('We can push rank %d to %d \n',i,r);
  %u=u(:,1:r);
  %v=v(:,1:r)*diag(s(1:r));
  v=v(:,1:r);
  v=reshape(v,[n2,ry3,r]); 
  v=permute(v,[1,3,2]);
  y{i+1}=v;
  %u=reshape(u,[ry1,n1,r]);
  %y{i}=permute(u,[2,1,3]);
  %v=reshape(v,[n2,ry3,r]);
  %y{i+1}=permute(v,[1,3,2]);
  %Finally, need to recompute ph{i}
  %will have size ry
  %A2(n2,m2,ra2,ra3)*X2(m1,rx2,rx3)*Y2(n2,ry2,ry3)*ph{i+2}(rx3,ra3,ry3)
  %ph_save(n2,ra2,rx2,ry3)*Y(n2,r,ry3]
  %Result is rx2,ra2,r, con
  ph_save=reshape(ph_save,[n2,ra2,rx2,ry3]);
  v=permute(v,[1,3,2]); v=reshape(v,[n2*ry3,r]);
  ph_save=permute(ph_save,[1,4,3,2]);
  ph_save=reshape(ph_save,[n2*ry3,rx2*ra2]);
  ph_save=ph_save'*v;
  ph{i+1}=reshape(ph_save,[rx2,ra2,r]);
 end
 
 
end 
 
 
 %And finally again recompute the first core
%A(n1,m1,ra1)*X(m1,rx1)*A(n2,m2,ra1,ra2)*X(m2,rx1,rx2)*ph(rx2,ra2,ry2)
%Output will be B(n1,n2,ry2), will separate n1; n2ry2


%Decipher sizes
core1=a{1}; n1=size(core1,1); m1=size(core1,2); ra1=size(core1,3);
core2=a{2}; n2=size(core2,2); m2=size(core2,2); ra2=size(core2,4);
x1=x{1}; rx1=size(x1,2); 
x2=x{2}; rx2=size(x2,3);
phi2=ph{3};
ry2=size(phi2,3);
%X(m2,rx1,rx2)*ph(rx2,ra2,ry2);
phi2=reshape(x2,[m2*rx1,rx2])*reshape(phi2,[rx2,ra2*ry2]);
phi2=reshape(phi2,[m2,rx1,ra2,ry2]);
phi2=permute(phi2,[1,3,2,4]);
core2=permute(core2,[1,3,2,4]);
phi2=reshape(core2,[n2*ra1,m2*ra2])*reshape(phi2,[m2*ra2,rx1*ry2]);
phi2=reshape(phi2,[n2,ra1,rx1,ry2]);
phi2=permute(phi2,[2,3,1,4]);
phi2=reshape(phi2,[ra1*rx1,n2*ry2]);
core1=permute(core1,[1,3,2]);
x1=reshape(core1,[n1*ra1,m1])*x1;
mat=reshape(x1,[n1,ra1*rx1])*phi2;
mat=reshape(mat,[n1,n2*ry2]);


[u,s,v]=svd(mat);
s=diag(s);
rold=size(y{d},2);
r=my_chop2(s,eps*norm(s));
r=max([r,rmin,rold+radd]);
r=min(r,size(s,1));
if ( verb )
fprintf('We can push rank %d to %d \n',1,r);
end
u=u(:,1:r);
v=v(:,1:r);
v=v*diag(s(1:r));
y{1}=u;
v=reshape(v,[n2,ry2,r]);
v=permute(v,[1,3,2]);
y{2}=v;
%v=v*diag(s(1:r));
%y{1}=u;
%y{2}=v; %v is n2*ry2*ry1
%v=reshape(v,[n2,ry2,r]);
%v=permute(v,[1,3,2]); v=reshape(v,[n2*r,ry2]); [v,dmp]=qr(v,0);
%v=reshape(v,[n2,r,ry2]);
%y{2}=v;
 return
end