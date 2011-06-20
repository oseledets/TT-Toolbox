%% HENON-HEILES POTENTIAL
a0=-7; b0=7;
q=1;
ev=1;
f=10;
for d=8:8
%a=tt_lapl(d);
lp=tt_lapl(d);
 n=2^d; h=(b0-a0)/(n+1);
x1=a0+(1:n)*h;
x=vec_to_nd(x1,1e-8,d);
%lp=ttm_add(lp,tt_vec_to_diag(x));

a=tt_qlaplace_dd(d*ones(1,f)); a=tt_matrix(a); a0=a; a=[];
a=a*0.5/h^2;
%Add Harmonic oscillator part;
xx=cell(f,1);
e=tt_ones(d,2); e=tt_tensor(e); x=tt_tensor(x);
%% The additional potential goes here
% lam*sum_{k=1}^{f-1} q^2_k * q_{k+1} - 1/3*q_{k+1}^3
 lam=0.111803;
 %lam=0;
for i=1:f
   xx{i}=x;
   for j=1:i-1       
       xx{i}=kron(e,xx{i});
   end
   for j=i+1:f
     xx{i}=kron(xx{i},e);
   end
end
ho=xx{1}.^2;
for i=2:f
    ho=ho+xx{i}.^2;
   ho=round(ho,1e-10);
end
ho=ho*0.5;
a=a+diag(ho); a=round(a,1e-10);
for i=1:f-1
    a=a+diag(lam*(xx{i}.*xx{i}.*xx{i+1}-1.0/3.0*xx{i+1}.^3));
   a=round(a,1e-10);
end
 %% Solution step
 tic;
 K=2;
%v=dmrg_eigb(a,K,1e-3,[],[],7);
%v1=a*v; ev1=dot(v1,v); norm(v1-ev1*v)

eps_loc=1e-4;
er1=2*eps_loc;
a=a0+a;
%while(er1>eps_loc)
%return
v=dmrg_eigb(a,5,1e-7,[],[],6);%v=v/norm(v);
%v1=a*v; ev1=dot(v1,v); 
%  er1=norm(v1-ev1*v);
%end
ev(q)=ev1; 
q=q+1;
toc;
end
return



