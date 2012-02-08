function [y,swp]=tt_mvk2(a,x,eps, max_r, max_swp)
%Matrix-by-vector product by DMRG/Krylov method
%   [A]=TT_MVK2(A,X,EPS, max_r, max_swp) Matrix-by-vector product by
%   DMRG/Krylov method. Try to use tt_mvk3/mvk3 methods: this one is likely
%   to be outdated and with not correct parameter selection
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
exists_max_r=1;
if ((nargin<4)||(isempty(max_r)))
    exists_max_r=0;
end;

nswp=15;
if ((nargin>=5)&&(~isempty(max_swp)))
    nswp = max_swp;
end;

kickrank = 5;
dropsweeps = 1; % garbage - for quasi-wedderburn
ddpow = 0.1; % stepsize for d-power in truncations
ddrank = 1; % stepsize for additional rank
d_pow_check = 0; % d-power for checking the convergence
bot_conv = 0.1; % bottom convergence factor - if better, we can decrease dpow and drank
top_conv = 0.99; % top convergence factor - if worse, we have to increase dpow and drank
verb = 2; % 0 - silent, 1 - sweep information, 2 - block information

d=size(a,1);

%Init solution
% a_coarse = tt_mat_compr(a, eps, 1);
% x_coarse = tt_compr2(x, eps, round(max(tt_ranks(x))/5));
% y = tt_mv(a_coarse, x_coarse);
% y=tt_mv(a_coarse,x);
% y = x;
%y=tt_mv(a,x);
y=tt_random(tt_size(a),d,2);
ph=cell(d,1);
%ph will be rx * ra * ry 
swp=1;

log_norms_ph = zeros(d,1);
log_norms_y = zeros(d,1);
last_sweep = false;

dy_old = ones(d-1,1);
dy = zeros(d-1,1);
% artificial rank additions
drank = ones(d-1,1);
% d-power for stronger compression eps./(d.^dpows)
dpows = ones(d-1,1);

while ( swp < nswp) 
 %power up by addition random rank-two term (it is magic)
%  y0=tt_random(tt_size(a),d,2);
%  y=tt_add(y,y0);
 
 for q=1:d
     log_norms_y(q)=log(norm(reshape(y{q}, size(y{q},1)*size(y{q},2), size(y{q},3)), 'fro')+1e-308);
     y{q}=y{q}./exp(log_norms_y(q));
 end;
 
   % fprintf('swp=%d \n',swp);
   err_max = 0;
    %Back reorthogonalization
  %Sweep starts from right-to-left orthogonalization of x
  %Compute right-to-left qr and phi matrices
  mat=y{d};
  [q,rv]=qr(mat,0);
  log_norms_y(d)=log(norm(q,'fro')+1e-308);
  y{d}=q;
  y{d}=y{d}./exp(log_norms_y(d));
  %A(n,m,ra))*x(m,rx)*y(n,ry) -> ph(rx,ra,ry)
  core=a{d}; n=size(core,1); m=size(core,2); ra=size(core,3);
  mat_x=x{d}; rx=size(mat_x,2);
  mat_y=y{d}; ry=size(mat_y,2);
  core=permute(core,[1,3,2]); 
  core=reshape(core,[n*ra,m]);
  core=core*mat_x; %is n*ra*rx1
  core=reshape(core,[n,ra,rx]); %n*ra*rxm
  core=permute(core,[2,3,1]); %ra*rxm*n
  core=reshape(core,[ra*rx,n]);
  ph{d}=core*mat_y; % ra*rx*ry
  log_norms_ph(d) = log(norm(ph{d}, 'fro')+1e-308); % preserve norms of phi to avoid Infs
  ph{d}=ph{d}./exp(log_norms_ph(d)); % keep only normalized phi
  log_norms_ph(d) = log_norms_ph(d) + log_norms_y(d);
  ph{d}=reshape(ph{d},[ra,rx,ry]);
  ph{d}=permute(ph{d},[2,1,3]);

%Back reorthogonalization
for i=(d-1):-1:2
    core=y{i};
    core=ten_conv(core,3,rv');
    ncur=size(core,1);
    r2=size(core,2);
    r3=size(core,3);
    core=permute(core,[1,3,2]);
    core=reshape(core,[ncur*r3,r2]);
    [y{i},rv]=qr(core,0);
    log_norms_y(i)=log(norm(y{i}, 'fro')+1e-308);
    y{i}=y{i}./exp(log_norms_y(i));    
    rnew=min(r2,ncur*r3);
    y{i}=reshape(y{i},[ncur,r3,rnew]);
    y{i}=permute(y{i},[1,3,2]);
    %A(n,m,ra1,ra2)*x(m,rx1,rx2)*y(n,ry1,ry2)*ph(rx2,ra2,ry2)->
    %phx(rx1,ra1,ry1)
    %Decipher sizes
    a1=a{i}; n=size(a1,1); m=size(a1,2); ra1=size(a1,3); ra2=size(a1,4);
    x1=x{i}; 
    rx1=size(x1,2); rx2=size(x1,3); 
    y1=y{i}; 
    ry1=size(y1,2); ry2=size(y1,3);
    phi=ph{i+1};
    x1=reshape(x1,[m*rx1,rx2]);
    phi=reshape(phi,[rx2,ra2*ry2]);
    phi=x1*phi; %ph is m*rx1*ra2*ry2
    phi=reshape(phi,[m,rx1,ra2,ry2]);
    %Convolve over m,ra2
    phi=permute(phi,[1,3,2,4]);
    phi=reshape(phi,[m*ra2,rx1*ry2]);
    a1=permute(a1,[1,3,2,4]);
    a1=reshape(a1,[n*ra1,m*ra2]);
    phi=a1*phi; %ph is n*ra1*rx1*ry2
    %Convolve over n*ry2
    phi=reshape(phi,[n,ra1,rx1,ry2]);
    phi=permute(phi,[1,4,2,3]);
    phi=reshape(phi,[n*ry2,ra1*rx1]);
    y1=permute(y1,[1,3,2]);
    y1=reshape(y1,[n*ry2,ry1]);
    phi=y1'*phi;
    log_norms_ph(i)=log(norm(phi, 'fro')+1e-308);
    phi = phi./exp(log_norms_ph(i)); % keep normalized phi
    log_norms_ph(i)=log_norms_ph(i)+log_norms_ph(i+1)+log_norms_y(i); % don't forget norm of previous phi
    %ph is ry1*ra1*rx1
    phi=reshape(phi,[ry1,ra1,rx1]);
    phi=permute(phi,[3,2,1]);
    ph{i}=phi;
%     norm_phi=norm(reshape(phi,ry1*ra1*rx1,1), 'fro')
end


%First, the first mode
%A(n1,m1,ra1)*X(m1,rx1)*A(n2,m2,ra1,ra2)*X(m2,rx1,rx2)*ph(rx2,ra2,ry2)
%Output will be B(n1,n2,ry2), will separate n1; n2ry2
%Decipher sizes


core1=a{1}; n1=size(core1,1); m1=size(core1,2); ra1=size(core1,3);
core2=a{2}; n2=size(core2,1); m2=size(core2,2); ra2=size(core2,4);
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

%Here is multiplication by the first core (to be substituted by a 
%sparse matrix-by-vector product)

% fprintf('norm_phi2=%g\n', norm(phi2));
core1=permute(core1,[1,3,2]);
x1=reshape(core1,[n1*ra1,m1])*x1;
% if ((max(max(isinf(x1)))==1))
%     fprintf('inf in x1\n');
% end;
% if ((max(max(isinf(phi2)))==1)||(max(max(isnan(phi2)))==1))
%     fprintf('inf in phi2\n');
% end;
mat=reshape(x1,[n1,ra1*rx1])*phi2;
mat=reshape(mat,[n1,n2*ry2]);

% if ((max(max(isinf(mat)))==1)||(max(max(isnan(mat)))==1))
%     fprintf('inf in mat\n');
% end;

[u,s,v]=svd(mat,0);
s=diag(s);
nrm=norm(s);
rold=size(y{1},2);
r=my_chop2(s,eps/sqrt(d)*nrm);
% r=max([rmin,r,rold+radd]);
% r=min(r,size(s,1));
if (exists_max_r) r = min(r, max_r); end;
if ( verb>1 ) 
 fprintf('We can push rank %d to %d \n',1,r);
end
u=u(:,1:r);
v=v(:,1:r)*diag(s(1:r));
% Power up by kickrank random addition
u = [u, randn(size(u,1),kickrank)];
v = [v, zeros(size(v,1),kickrank)];
[u,rv]=qr(u,0);
r = size(u,2);
v = v*(rv.');
y{1}=u;
log_norms_y(1)=log(norm(u, 'fro')+1e-308);
y{1}=y{1}./exp(log_norms_y(1));
log_norms_y(2)=log(norm(v, 'fro')+1e-308);
y{2}=reshape(v,[n2,ry2,r]);
y{2}=permute(y{2},[1,3,2]);
y{2}=y{2}./exp(log_norms_y(2));
log_norms_y(2)=log_norms_y(2)+log_norms_ph(3); % v*s' was generated using ph{3}
%Recalculate ph{1};
%ph{1}=A(n,m,ra)*X(m,rx)*Y(n,ry) -> (ry,ra,rx)
%Decipher sizes
core=a{1}; n=size(core,1); m=size(core,2); ra=size(core,3);
x1=x{1}; rx=size(x1,2);
y1=y{1}; ry=size(y1,2);
core=permute(core,[1,3,2]); core=reshape(core,[n*ra,m]);
phi1=core*x1; %phi is n*ra*rx
phi1=y1'*reshape(phi1,[n,ra*rx]);
log_norms_ph(1)=log(norm(phi1, 'fro')+1e-308);
phi1 = phi1./exp(log_norms_ph(1));
log_norms_ph(1) = log_norms_ph(1) + log_norms_y(1); % phi1=y1*core*x1
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
  %Check with old approximation
  y1=y{i};
  y2=y{i+1};
  ry2=size(y1,3);
  y1=permute(y1,[2,1,3]); 
  y2=permute(y2,[2,1,3]);
  y1=reshape(y1,[ry1*n1,ry2]);
  y2=reshape(y2,[ry2,n2*ry3]);
  y2 = y2*exp(log_norms_y(i)+log_norms_y(i+1)-log_norms_ph(i-1)-log_norms_ph(i+2)); % recover the true norm of y2
  app=y1*y2;
  dy(i)=norm(mat-app,'fro')/norm(mat,'fro');
  if (swp==1)
      dy_old(i)=dy(i);
  end;
  
  % The new core does not converge - increase rank
  if (dy(i)/dy_old(i)>top_conv)&&(dy(i)>eps/(d^d_pow_check))
      drank(i)=drank(i)+ddrank;
      dpows(i)=dpows(i)+ddpow;
  end;
  % The new core converges well - try to decrease rank
  if (dy(i)/dy_old(i)<bot_conv)||(dy(i)<eps/(d^d_pow_check))
      drank(i)=max(drank(i)-ddrank, 1);
      dpows(i)=max(dpows(i)-ddpow, 1);
  end;
  % perform simple compression for the last sweep
  if (last_sweep)
      dpows(i)=0.5;
  end;
  
%   if ( er0 > eps )
% %     fprintf('rank %d\t does not converge,er0=%3.2e\t\n',i,er0);
%     not_converged=true;
%   end
  %fprintf('er0=%3.2e\n',er0);
  if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
      [u,s,v]=svd(mat-app,'econ');
  else
      [u,s,v]=svd(mat,'econ');
  end;
   s=diag(s);
  
  %fprintf('norm=%18f \n',norm(s));
  %fprintf('tensor norm=%3.2e \n',norm(s));

  r=my_chop2(s,eps/(d^dpows(i))*norm(mat, 'fro'));
  if (~last_sweep)
      r = r+drank(i); % we want even larger ranks
  end;
  r = min(r, max(size(s))); % but not too large
  if (exists_max_r) r = min(r, max_r); end;  
  %if ( 
  if ( verb>1 )
  fprintf('We can push rank %d to %d \n',i,r);
  end
  u=u(:,1:r);
  v=conj(v(:,1:r))*diag(s(1:r));
  if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
      u = [y1, u];
      v = [y2.', v];
      [u,rv]=qr(u,0);
      v = v*(rv.');
      r = size(u,2);
  else
      if (~last_sweep)
          % kick
          u = [u, randn(size(u,1),kickrank)];
          v = [v, zeros(size(v,1),kickrank)];
          [u,rv]=qr(u,0);
          r = size(u,2);
          v = v*(rv.');
      end;
  end;
  log_norms_y(i)=log(norm(u, 'fro')+1e-308);  
  log_norms_y(i+1)=log(norm(v, 'fro')+1e-308);  
  u=reshape(u,[ry1,n1,r]);
  y{i}=permute(u,[2,1,3]);
  y{i}=y{i}./exp(log_norms_y(i));
  v=reshape(v,[n2,ry3,r]);
  y{i+1}=permute(v,[1,3,2]);
  y{i+1}=y{i+1}./exp(log_norms_y(i+1));
  log_norms_y(i+1)=log_norms_y(i+1)+log_norms_ph(i-1)+log_norms_ph(i+2); % norm(mat) ~ norm(ph{i-1})*norm(ph{i+2})
  %Finally, need to recompute ph{i}
  %will have size ry
  %A(n1,m1,ra1,ra2)*X(m1,rx1,rx2)*Y(n1,ry1,ry2)*ph{i-1}(ry1,ra1,rx1)
  %ph_save(rx2,ry1,ra2,n1)*Y(n1,ry1,ry2) 
  ph_save=reshape(ph_save,[rx2,ry1,ra2,n1]);
  ph_save=permute(ph_save,[4,2,3,1]);
  ph_save=reshape(ph_save,[n1*ry1,ra2*rx2]);  %HERE WE CAN HAVE A FAULT
  ph_save=reshape(y{i},[n1*ry1,r])'*ph_save;
  log_norms_ph(i)=log(norm(ph_save, 'fro')+1e-308);
  ph_save = ph_save./exp(log_norms_ph(i));
  log_norms_ph(i)=log_norms_ph(i)+log_norms_y(i)+log_norms_ph(i-1);
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
 r=my_chop2(s,eps/sqrt(d)*norm(s));
 rold=size(y{d},2);
% r=max([r,rmin,rold]);
% r=min(r,size(s,1));
 if (exists_max_r) r = min(r, max_r); end;
 %fprintf('We can push rank %d to %d \n',d-1,r);
 u=u(:,1:r);
 v=v(:,1:r)*diag(s(1:r)); 
 % kick
 u = [u, randn(size(u,1),kickrank)];
 v = [v, zeros(size(v,1),kickrank)];
 [u,rv]=qr(u,0);
 r = size(u,2);
 v = v*(rv.');
 %v=v;
 log_norms_y(d-1)=log(norm(u, 'fro')+1e-308);
 log_norms_y(d)=log(norm(v, 'fro')+1e-308);
 u=reshape(u,[ry1,n1,r]);
 u=permute(u,[2,1,3]);
 y{d-1}=u;
 y{d-1}=y{d-1}./exp(log_norms_y(d-1));
 y{d}=v;
 y{d}=y{d}./exp(log_norms_y(d));
 log_norms_y(d-1)=log_norms_y(d-1)+log_norms_ph(d-2);
 %Recalculate ph{d}
 %A2(n2,m2,ra2)*X(m2,rx2)*Y(n2,ry2} %ph_save is n2*ra2*x2, result
 %rx2*ra2*ry2
 ph_save=reshape(ph_save,[n2,ra2,rx2]);
 ph_save=permute(ph_save,[1,3,2]); ph_save=reshape(ph_save,[n2,rx2*ra2]);
 ph_save=ph_save'*v;
 log_norms_ph(d)=log(norm(ph_save, 'fro')+1e-308);
 ph_save = ph_save./exp(log_norms_ph(d));
 ph_save=reshape(ph_save,[rx2,ra2,r]);
 ph{d}=ph_save;

%  norms_phi = exp(log_norms_ph')
%  norms_y = exp(log_norms_y')
 
 norm1_y = exp(sum(log_norms_y)/d);
 for i=1:d
     y{i}=y{i}.*norm1_y;
 end;
 
 if (verb>0)
     fprintf('tt_mvk2: sweep %d, err_max = %3.3e\n', swp, max(dy));
 end; 
 if (last_sweep)
     swp=swp+1;
     break;
 end;
 if (max(dy)<=eps/(d^d_pow_check))
     last_sweep=true;
 end;
 
 dy_old = dy;
  swp=swp+1;
  if (swp==nswp-1)
      last_sweep=true;
  end;
end
if ( swp == nswp )&&(max(dy)>eps/(d^d_pow_check)) 
  fprintf('tt_mvk2 warning: error is not fixed for maximal number of sweeps %d, err_max: %3.3e\n', swp, err_max); 
end

return
end
