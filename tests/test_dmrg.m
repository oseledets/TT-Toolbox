d=9;
d_phys=3;
% q=1000;
eps=1e-6;
maxit=20;

% [x,y]=meshgrid((1:1:2^d)/(2^d+1));
% q=-exp(-x.^2-y.^2)*500;
% q=-0*ones(2^d,2^d); % sqrt((x-0.5).^2+(y-0.5).^2);
% q=-1000*ones(1024,1024);
% q=-400*(x+y);
% q(1:512,1:512)=-400;
% q(:,1:512)=0;
x=(1:1:2^d)'/(2^d+1);
tt_x=full_to_qtt(x, eps);
tt_y = tt_kron_dc(tt_x, tt_ones(d,2),0);
tt_y = tt_kron_dc(tt_ones(d,2), tt_y, 0);
tt_z = tt_kron_dc(tt_x, tt_ones(d*2,2),0);
tt_x = tt_kron_dc(tt_ones(d*2,2), tt_x, 0);
% q=exp(-(x-0.5).^2);
tt_q = tt_add(tt_x, tt_y);
tt_q = tt_scal(tt_q, -400);
% tt_q=full_to_qtt(q, eps, d*d_phys);
% tt_q=full_to_qtt(-400*(1:1:2^d)'/(2^d+1), eps);
% tt_q=tt_q0;
% for i=2:d_phys
%     tt_q=tt_kron_dc(tt_q, tt_q0, 0);
% end;
% tt_q=tt_scal(tt_q, 100);

% a = 100*sqrt((x-0.5).^2+(y-0.5).^2)+1;
% a=100*exp(x.^2.*y)+1;
% tt_a=full_to_qtt(a, eps, d*d_phys);
% tt_a=tt_ones(d_phys, 2^d);
% A=tt_Fd_mtx(d_phys, tt_a, 0, 0, eps);
% tt_q=tt_ones(d*d_phys,2);
A=tt_qlaplace_dd(d*ones(1,d_phys));
A=tt_scal(A, (2^d+1)^2);
A=ttm_add(A, tt_vec_to_diag(tt_q));

% A2=ttm_add(tt_kron3(A, tt_eye(2,d*d_phys)), tt_kron3(tt_eye(2,d*d_phys), A));

A=tt_mat_compr(A, eps);

% par_values = cell(2*d_phys,1);
% for i=1:2:d_phys*2
%     par_values{i}=(0:1:2^d-1)'/(2^d-1);
%     par_values{i+1}=(0:1:2^d-1)'/(2^d-1);
% end;
% par_values{1}=0.5+(0:1:2^d-1)'*0.5/(2^d-1);
% A = tt_bound_param_Fd_mtx(tt_a, par_values, eps);
% A=tt_to_qtt3(A, 1, eps);

rhs=tt_ones(d_phys*d, 2);
% rhs=tt_mat_to_vec(tt_scal(tt_eye(2,d*d_phys), 2));
x=tt_random(2, d_phys*d, 2);
% x=tt_scal(rhs,0.5);

results=zeros(maxit, 4);

for i=1:maxit
    tic;
    x=dmrg_solve2(A, rhs, x, eps, eps, 100, 1);
    cur_time=toc
    max_xrank=max(tt_ranks(x))

    resid=tt_dist2(tt_mv(A,x),rhs)/sqrt(tt_dot(rhs,rhs))

    
    results(i,1)=i;
    results(i,2)=resid;
    results(i,3)=max_xrank;
    results(i,4)=cur_time;
    pause(2);
end;

fulltimes=cumsum(results(:,4));
for j=1:i
    fprintf('%d\t%3.3e\t%d\t%g\t%g\n', results(j,1), results(j,2), results(j,3), results(j,4), fulltimes(j));
end;

% for d=9:10
%     lp1=tt_qlaplace_dd(d);
%     lp1=tt_scal(lp1, (2^d+1)^2);
%     I=tt_eye(2,d);
%     
%     A=ttm_add(tt_kron3(lp1, I), tt_kron3(I, lp1));
% %     A=ttm_add(tt_kron_dc(lp1,I,1), tt_kron_dc(I, lp1, 1));
%     % A=ttm_add(A, tt_kron3(tt_kron3(I,I), lp1));
%     % A = tt_qlaplace_dd([d,d,d]);
%     % A=tt_scal(A, (2^d+1)^2);
%     A=ttm_add(A, tt_scal(tt_eye(4,d), -q));
%     A=tt_mat_compr(A, eps);
%     
%     rhs=tt_ones(d, 4);
%     
%     if (d==9) x0=rhs; end;
%     x=dmrg_solve(A, rhs, x0, eps, eps, 100, 3)
%     
%     tt_dist2(tt_mv(A,x),rhs)/sqrt(tt_dot(rhs,rhs))
%     
%     prol0 = [1 0; 0.5 0.5; 0 1; 0 0.5]; % 4x2
%     h0 = [0 0; 0 0; 0 0; 0.5 0];
%     I = eye(2);
%     I12 = [0 1; 0 0];
%     I21 = [0 0; 1 0];
%     prol=cell(d,1);
%     prol{d} = zeros(4,4,2);
%     prol{d}(:,:,1)=kron(I,I); prol{d}(:,:,2)=kron(I12, I12);
%     prol{1}=zeros(16,4,2);
%     prol{1}(:,:,1)=kron(prol0,prol0); prol{1}(:,:,2)=kron(h0,h0);
%     for i=2:d-1
%         prol{i}=zeros(4,4,2,2);
%         prol{i}(:,:,1,1)=kron(I,I);
%         prol{i}(:,:,2,1)=kron(I12,I12); prol{i}(:,:,2,2)=kron(I21,I21);
%     end;
%     
%     x2=tt_mv(prol,x);
%     rv=1; rnew=1;
%     for i=d:-1:2
%         core=x2{i};
%         n1=size(core,1); r1=size(core,2); r2=size(core,3);
%         core=reshape(core, n1*r1, r2);
%         core=core*rv.'; % size n1*r1, r2new
%         core=permute(reshape(core, n1, r1, rnew), [1 3 2]);
%         r2=rnew;
%         core=reshape(core, n1*r2, r1);
%         [qnew,rv]=qr(core,0); % size n1*r2, rnew - rnew,r1
%         rnew=min(n1*r2, r1);
%         x2{i}=permute(reshape(qnew, n1, r2, rnew), [1 3 2]);
%     end;
%     core=x2{1};
%     n1=size(core,1); r2=size(core,2);
%     core=core*rv.'; % size n1, r2new
%     % n1=16 - 4(x1)*4(x2). We need (2(x1)*2(x2))*(2(x1)*2(x2))
%     core=reshape(core, 2,  2,  2,  2, rnew);
%     %                 x10 x11 x20 x21
%     core=permute(core, [1 3 2 4 5]); % x10*x20, x11*x21*rnew
%     core=reshape(core, 4, 4*rnew);
%     [u,s,v]=svd(core,'econ');
%     r=size(u,2);
%     core0=reshape(u*s, 4, r);
%     core1=permute(reshape(v', r, 4, rnew), [2 1 3]);
%     
%     x0=cell(d+1,1);
%     for i=2:d
%         x0{i+1}=x2{i};
%     end;
%     x0{2}=core1;
%     x0{1}=core0;
%     
% %     x2=cell(d+1,1);
% %     for i=1:d
% %         x2{i+1}=x{i};
% %     end;
% %     x2{1}=ones(4,1);
% %     x2{2}=reshape(x2{2}, 4, 1, size(x2{2},2));
%     
% %     x0=x2;
% end;