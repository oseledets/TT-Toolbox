function [x, somedata]=dmrg_solve2(A, y, eps,varargin)
%Solution of linear systems in TT-format via DMRG iteration
%   [X,SWEEPS]=DMRG_SOLVE2(A,Y,EPS,OPTIONS) Attempts to solve the linear
%   system A*X = Y with accuracy EPS using the two-sided DMRG iteration.
%   Matrix A has to be given in the TT-format, right-hand side Y should be
%   given in the TT-format also. Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following) 
%   The list of option names and default values are:
%       o x0 - initial approximation [random rank-2 tensor] 
%       o P - preconditioner  [I]
%       o nswp - maximal number of DMRG sweeps [10]
%       o rmax - maximal TT-rank of the solution [1000]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o max_full_size - maximal size of the local matrix to full solver 
%       [2500]
%       o local_prec: Local preconditioner, 'als' - ALS-Richardson
%       iteration, 'selfprec' (Saad selfpreconditioner) ['als']
%       o prec_compr - compression for local precs [1e-3]
%       o prec_tol - tolerance for local precs [1e-1]
%       o prec_iters - number of local iterations [15]
%       o use_self_prec - Use self_prec [ true | {false} ]
%       o gmres_iters - number of local gmres restarts [2]
%       o nrestart - dimension of local gmres [40]
%       Example:
%           d=8; f=8; 
%           mat=tt_qlaplace_dd(d*ones(1,f)); %Laplace in the QTT-format
%           rhs=tt_ones(2,d*f); Right-hand side of all ones
%           sol=dmrg_solve2(mat,rhs,1e-6);
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


% Inner parameters
max_full_size=2500;
prec_compr=1e-3;
prec_tol=1e-1;
prec_iters=10;

dropsweeps=1;
ddpow = 0.1; % stepsize for d-power in truncations
min_dpow = 1; % Minimal d-power for truncation
ddrank = 1; % stepsize for additional rank
min_drank = 1; % Minimal additional rank
d_pow_check = 0; % d-power for checking the convergence
bot_conv = 0.1; % bottom convergence factor - if better, we can decrease dpow and drank
top_conv = 0.99; % top convergence factor - if worse, we have to increase dpow and drank

bs_treshold = 0.0001*0; % Treshold from the previous residual to consider a local system as "bad"
trunc_to_true = 2; % Truncation error to true residual treshold

use_self_prec=false;
nswp=10;
nrestart=40;
gmres_iters=2;
% local_prec = 'als';
local_prec = 'selfprec';
local_format = 'full';
% local_format = 'tt';
rmax=1000;
tol=eps;
verb=1;
kickrank = 2;
x0=[];
P=[];
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
            rmax=lower(varargin{i+1});
        case 'x0'
            x0=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'p'
            P=varargin{i+1};
        case 'tol'
            tol=varargin{i+1};
        case 'local_prec'
            local_prec=varargin{i+1};
        case 'nrestart'
            nrestart=varargin{i+1};
        case 'gmres_iters'
            gmres_iters=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case  'max_full_size'
            max_full_size=varargin{i+1};
        case 'prec_compr'
            prec_compr=varargin{i+1};
        case 'prec_tol'
            prec_tol=varargin{i+1};
        case 'prec_iters'
            prec_iters=varargin{i+1};
        case 'use_self_prec'
            use_self_prec=varargin{i+1};
        case 'ddpow'
            ddpow=varargin{i+1};
        case 'ddrank'
            ddrank=varargin{i+1};
        case 'd_pow_check'
            d_pow_check=varargin{i+1};
        case 'bot_conv'
            bot_conv=varargin{i+1};
        case 'top_conv'
            top_conv=varargin{i+1};
        case 'min_dpow'
            min_dpow=varargin{i+1};
        case 'min_drank'
            min_drank=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

input_is_tt_tensor = 0;

if ( isa(y,'tt_tensor') )
  y=core(y);
  input_is_tt_tensor = 1;
end

if (isa(y, 'tt_tensor'))
    y=core(y);
    input_is_tt_tensor = 1;
end;
if ( isa(A,'tt_matrix') )
%   ttA=A.tt;
%   dA=ttA.d;
  A=core(A);
  input_is_tt_tensor = 1;
  %if (isempty(x0))
  %    x0=tt_random(tt_size(y), A.tt.d, 2);
  %end;
% else
%    dA=numel(A);
   % if (isempty(x0))
   %     x0=tt_random(tt_size(y), max(size(A)), 2);
   % end;
end
%   x0=tt_random(tt_size(y),dA,2);
if (isempty(x0))
    x0=tt_random(tt_size(y), max(size(y)), 2);
end;
if ( isa(P,'tt_matrix') )
  P=core(P);
  input_is_tt_tensor = 1;
end;
if ( isa(x0,'tt_tensor') )
  x0=core(x0);
  input_is_tt_tensor = 1;
end

%nrmF=sqrt(tt_dot(y,y));

d=size(A,1);

if ( isempty(P) )
   P = core(tt_eye(tt_size(y), d));
end
x=x0;

x{1}=reshape(x{1}, size(x{1},1), 1, size(x{1},2));
y{1}=reshape(y{1}, size(y{1},1), 1, size(y{1},2));
A{1}=reshape(A{1}, size(A{1},1),size(A{1},2), 1, size(A{1},3)); %Bydlocode (@)
P{1}=reshape(P{1}, size(P{1},1),size(P{1},2), 1, size(P{1},3));

phA = cell(d,1);
phy = cell(d,1);
dx_old = ones(d,1);
dx = zeros(d,1);
% artificial rank additions
drank = ones(d,1)*min_drank;
% d-power for stronger compression eps./(d.^dpows)
dpows = ones(d,1)*min_dpow;

%  chkvec = tt_random(tt_size(y), max(size(y)), kickrank);
%  chkvec{1}=reshape(chkvec{1}, size(chkvec{1},1), 1, size(chkvec{1},2));
%  phAchk = cell(d,1);
%  phychk = cell(d,1);


somedata = cell(4,1); % swp, conds
somedata{2} = zeros(d, nswp);

sol_hist = cell(3,1);
sol_hist{1}=x;
max_res_old = 0;
last_sweep = false;
for swp=1:nswp
%     z = x;
    % 1-to-d orthogonalization
    rvx = 1; rnewx=1; phAold=1; phyold=1;

    for i=1:d-1
        cre = x{i};
        n1 = size(cre,1); rx1 = size(cre,2); rx2 = size(cre,3);
        cre = reshape(permute(cre, [2 1 3]), rx1, n1*rx2);
        cre = rvx*cre; % size rnew,n1,rx2
        rx1=rnewx;
        cre = reshape(cre, rx1*n1, rx2);
        [q,rvx]=qr(cre,0); % size rx1*n1,r2new - r2new,rx2
        rnewx = min(rx1*n1, rx2);
        x{i}=permute(reshape(q, rx1, n1, rnewx), [2 1 3]);

        % Now, update phi. phA=X' PA X, phY = X' PY
        a1 = A{i};
        n1=size(a1,1); m1=size(a1,2); ra1=size(a1,3); ra2=size(a1,4);
        p1 = P{i};
        k1=size(p1,1); rp1=size(p1,3); rp2=size(p1,4);
        y1 = y{i};
        ry1=size(y1,2); ry2=size(y1,3);
        x1 = x{i};

        rxm1=size(x1,2); rxm2=size(x1,3); rxn1=rxm1; rxn2=rxm2;
        phAold = reshape(phAold, rxn1*rp1*ra1, rxm1);
        x1 = reshape(permute(x1, [2 1 3]), rxm1, m1*rxm2);
        phAold=phAold*x1; % size rxn1*rp1*ra1*m1*rxm2
        phAold=reshape(phAold, rxn1, rp1, ra1, m1, rxm2);
        phAold=reshape(permute(phAold, [1 2 5 4 3]), rxn1*rp1*rxm2, m1*ra1);
        a1 = reshape(permute(a1, [2 3 1 4]), m1*ra1, n1*ra2);
        phAold=phAold*a1; % size rxn1*rp1*rxm2*n1*ra2
        phAold=reshape(phAold, rxn1, rp1, rxm2, n1, ra2);
        phAold=reshape(permute(phAold, [1 3 5 4 2]), rxn1*rxm2*ra2, n1*rp1);
        p1 = reshape(permute(p1, [2 3 1 4]), n1*rp1, k1*rp2);
        phAold=phAold*p1; % size rxn1*rxm2*ra2*k1*rp2
        phAold=reshape(phAold, rxn1, rxm2, ra2, k1, rp2);
        phAold=reshape(permute(phAold, [2 5 3 4 1]), rxm2*rp2*ra2, k1*rxn1);
        x1 = reshape(x{i}, k1*rxn1, rxn2);
        phAold=phAold*conj(x1); % size rxm2*rp2*ra2*rxn2 <--- complex conjugate!
        phAold = reshape(phAold, rxm2, rp2, ra2, rxn2);
        phAold = permute(phAold, [4 2 3 1]); % we need rxn2,rp2,ra2,rxm2
        phA{i}=phAold;

        phyold = reshape(phyold, rxn1*rp1, ry1);
        y1 = reshape(permute(y1, [2 1 3]), ry1, n1*ry2);
        phyold=phyold*y1; % size rxn1*rp1*n1*ry2
        phyold = reshape(phyold, rxn1, rp1, n1, ry2);
        phyold=reshape(permute(phyold, [1 4 3 2]), rxn1*ry2, n1*rp1);
        p1=reshape(permute(P{i}, [2 3 1 4]), n1*rp1, k1*rp2);
        phyold=phyold*p1; % size rxn1*ry2*k1*rp2
        phyold=reshape(phyold, rxn1, ry2, k1, rp2);
        phyold=reshape(permute(phyold, [4 2 3 1]), rp2*ry2, k1*rxn1);
        x1=reshape(x{i}, k1*rxn1, rxn2);
        phyold=phyold*conj(x1); % size rp2*ry2*rxn2 <--- complex conjugate!
        phyold=permute(reshape(phyold, rp2, ry2, rxn2), [3 1 2]);
        phy{i}=phyold;
    end;
    % convolve rv with the last cre
    cre = x{d};
    n1 = size(cre,1); rx1 = size(cre,2); rx2 = size(cre,3);
    cre = reshape(permute(cre, [2 1 3]), rx1, n1*rx2);
    cre = rvx*cre; % size rnew,n1,rx2
    x{d}=permute(reshape(cre, rnewx, n1, rx2), [2 1 3]);

    % Now, start the d-to-1 DMRG iteration
    max_res = 0;
    phAold=1; phyold=1;

    for i=d:-1:2
        a2=A{i}; a1=A{i-1}; ra1=size(a1,3); ra2=size(a1,4); ra3=size(a2,4);
        n1 = size(a1,1); m1=size(a1,2); n2=size(a2,1); m2=size(a2,2);
        p2=P{i}; p1=P{i-1}; rp1=size(p1,3); rp2=size(p1,4); rp3=size(p2,4);
        k1 = size(p1,1); k2=size(p2,1);

        y1=y{i-1}; y2=y{i}; ry1=size(y1,2); ry2=size(y1,3); ry3=size(y2,3);
        x1=x{i-1}; x2=x{i}; rx1=size(x1,2); rx2=size(x1,3); rx3=size(x2,3);

        % Compute RHS: phy{i-2}*P1*y1*y2*P2*phyold
        if (i>2)
            rhs1 = phy{i-2};
        else
            rhs1=1;
        end;
        rhs1 = reshape(rhs1, rx1*rp1, ry1);
        y1 = reshape(permute(y1, [2 1 3]), ry1, n1*ry2);
        rhs1 = rhs1*y1; % size rx1*rp1*n1*ry2
        rhs1 = reshape(rhs1, rx1, rp1, n1, ry2);
        rhs1=reshape(permute(rhs1, [1 4 3 2]), rx1*ry2, n1*rp1);
        p1 = reshape(permute(p1, [2 3 1 4]), n1*rp1, k1*rp2);
        rhs1=rhs1*p1; % size rx1*ry2*k1*rp2
        rhs1=reshape(rhs1, rx1, ry2, k1, rp2);
        rhs1=reshape(permute(rhs1, [1 3 4 2]), rx1*k1, rp2*ry2);

        y2=reshape(permute(y2, [2 1 3]), ry2*n2, ry3);
        phyold2 = reshape(phyold, rx3*rp3, ry3);
        rhs2 = y2*(phyold2.'); % size ry2*n2, rx3*rp3
        rhs2 = reshape(rhs2, ry2, n2, rx3, rp3);
        rhs2 = permute(rhs2, [1 3 2 4]);
        rhs2 = reshape(rhs2, ry2*rx3, n2*rp3);
        p2 = reshape(permute(p2, [2 4 1 3]), n2*rp3, k2*rp2);
        rhs2 = rhs2*p2; % size ry2*rx3, k2*rp2
        rhs2 = reshape(rhs2, ry2, rx3, k2, rp2);
        rhs2 = permute(rhs2, [4 1 3 2]);
        rhs2 = reshape(rhs2, rp2*ry2, k2*rx3);

        if (strcmp(local_format, 'full'))
            rhs = rhs1*rhs2;
            rhs = reshape(rhs, rx1*k1*k2*rx3, 1);
        else
            rhs = cell(2,1);
            rhs{1} = rhs1;
            rhs{2} = rhs2.';
        end;

        rxn1=rx1; rxn3=rx3;
        rxm1=rx1; rxm2=rx2; rxm3=rx3;
        if (i>2)
            B = phA{i-2};
        else
            B=1;
        end;
              
        if (i>2)
            B = phA{i-2};
        else
            B=1;
        end;        
        
        B = reshape(B, rxn1, rp1, ra1, rxm1);
        B = reshape(permute(B, [1 4 2 3]), rxn1*rxm1*rp1, ra1);
        a1 = reshape(permute(A{i-1}, [3 2 1 4]), ra1, m1*n1*ra2);
        B = B*a1; % size rxn1*rxm1*rp1*m1*n1*ra2
        B = reshape(B, rxn1,rxm1,rp1,m1,n1,ra2);
        B = permute(B, [1 2 4 6 3 5]);
        B = reshape(B, rxn1*rxm1*m1*ra2, rp1*n1);
        p1 = reshape(permute(P{i-1}, [3 2 1 4]), rp1*n1, k1*rp2);
        B = B*p1; % size rxn1*rxm1*m1*ra2*k1*rp2
        B = reshape(B, rxn1,rxm1,m1,ra2,k1,rp2);
        B = permute(B, [1 5 2 3 6 4]);
        B = reshape(B, rxn1*k1*rxm1*m1, rp2*ra2);
%         This is the first term of tensor-structured matrix B \otimes B2
%         Now, the second
        B2 = permute(phAold, [1 2 4 3]);
        B2 = reshape(B2, rxn3*rp3*rxm3, ra3);
        a2 = reshape(A{i}, n2*m2*ra2, ra3);
        B2 = B2*(a2.'); % size rxn3*rp3*rxm3*n2*m2*ra2
        B2 = reshape(B2, rxn3, rp3, rxm3, n2, m2, ra2);
        B2 = permute(B2, [1 3 5 6 4 2]);
        B2 = reshape(B2, rxn3*rxm3*m2*ra2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        B2 = B2*p2; % size rxn3*rxm3*m2*ra2*k2*rp2
        B2 = reshape(B2, rxn3, rxm3, m2, ra2, k2, rp2);
        B2 = permute(B2, [5 1 3 2 6 4]);
        B2 = reshape(B2, k2*rxn3*m2*rxm3, rp2*ra2);
        % Now, compress inner rank rp2*ra2 --- what is this ???
        %Modify it by random noise, since sometime MATLAB QR
        %fails
        %B2=B2+max(abs(B2(:)))*randn(size(B2))*1e-16; Kill all humans for
        %such code

%         [Q,R]=qr(B2,0);
%
%         rnew = min(k2*rxn3*m2*rxm3, rp2*ra2);
%         B2 = reshape(Q, k2*rxn3*m2*rxm3, rnew);
%         B = B*(R.'); % size rxn1*rxm1*k1*m1*rnew
%
%         B = reshape(B, rxn1*k1*rxm1*m1, rnew);
%         [U,S,V]=svd(B, 'econ');
%         S = diag(S);
%         rB = my_chop2(S, 1e-12*norm(S)); % We don't know the cond(B), so let's obtain almost exact compression
%         B = U(:,1:rB);
%         V = V(:,1:rB)*diag(S(1:rB)); % size rnew*rB
%         B2 = B2*conj(V); % size k2*rxn3*m2*rxm3*rB
        rB=rp2*ra2;
        MatVec='bfun2';
        if (((rxn1*k1*k2*rxn3<max_full_size))||(rB>max(rxn1*k1, rxn3*k2)))&&(strcmp(local_format, 'full'))
            MatVec='full';
            if (rxn1*k1*k2*rxn3>max_full_size)
                MatVec='half-full';
            end;
            B = B*(B2.'); % size rxn1*k1*rxm1*m1*k2*rxn3*m2*rxm3
            B = reshape(B, rxn1,k1,rxm1,m1,k2,rxn3,m2,rxm3);
            B = permute(B, [1 2 5 6 3 4 7 8]);
            B = reshape(B, rxn1*k1*k2*rxn3, rxm1*m1*m2*rxm3);
        else
            B1 = reshape(B, rxn1*k1, rxm1*m1, rB);
            B2 = reshape(B2, k2*rxn3, m2*rxm3, rB);
            B=cell(2,1);
            B{1}=B1;
            B{2}=B2;
        end;

        % Form previous solution
        x1 = reshape(permute(x{i-1}, [2 1 3]), rxm1*m1, rxm2);
        x2 = reshape(permute(x{i}, [2 1 3]), rxm2, m2*rxm3);
        if (strcmp(local_format, 'full'))
            sol_prev = x1*x2;
            sol_prev = reshape(sol_prev, rxm1*m1*m2*rxm3, 1);
        else
            sol_prev = cell(2,1);
            sol_prev{1} = x1;
            sol_prev{2} = x2.';
        end;

        real_tol = (tol/(d^dpows(i)))/trunc_to_true;

        if (strcmp(local_format, 'tt'))
            mv = @(vec,eps,mr)bfun3(B,vec,eps,mr);
        else
            if (strcmp(MatVec, 'bfun2'))
                mv=@(vec)bfun2(B, vec, rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
                mv_t=@(vec)bfun2_t(B, vec, rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
                %mv1=@(vec)bfun2(B, vec, rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)+tau(i)*vec;
            else
                mv = @(vec)(B*vec);
                mv_t = @(vec)(B'*vec);
                %mv1 = @(vec)(B*vec+tau(i)*vec);
            end;
        end;

        % Check the previous residual
        if (strcmp(local_format, 'tt'))
            res_prev = mv(sol_prev, [], []);
            normf = exp(0.5*tt_dot2(rhs, rhs));
            res_prev = tt_dist3(res_prev, rhs)/normf;
        else
            res_prev = mv(sol_prev);
            normf = norm(rhs);
            res_prev = norm(res_prev - rhs)/normf;
        end;

        % We will solve the system only if res_prev>0.1*max_res_prev
        if (~last_sweep)&&(res_prev>bs_treshold*max_res_old)
            if (strcmp(MatVec,'full'))
                
%                 ev=eig(B);
%                 ev1 = max(ev);
%                 evn = min(ev);
%                 [v1,ev1]=eigs(B'*B, [], 1, 'lr');
%                 [vn,evn]=eigs(B'*B, [], 1, 'sr');                
%                 somedata{2}(i, swp) = cond(B);
                somedata{3}(i, swp) = res_prev;
%                 fprintf('i=%d, cond(B): %g\n', i, somedata{2}(i, swp));
%                 keyboard;
                %             sol = pinv(B)*rhs;
                %             sol = (B'*B+tol^2*max(max(abs(B'*B)))*eye(size(B)))\(B'*rhs);
                sol = B \ rhs;
                %             sol = (B'*B)\(B'*rhs);
                res=B*sol;
                res_true = norm(res-rhs)/norm(rhs);
            else
                %Ax_{k+1}+tau*x_k=rhs+tau*x_k
                %(Ax_{k+1}+tau*I)x_{k+1}=rhs+tau*x_k
%                 [sol_new,flg] = gmres(mv1, rhs+tau(i)*sol_prev, nrestart, real_tol, 2, [], [], sol_prev);
%                 if( flg == 0)
%                     tau(i)=tau(i)/10;
%                 else
%                     tau(i)=tau(i)*4;
%                 end

%                    keyboard;
%                 [v1,ev1]=eigs(@(vec)(mv_t(mv(vec))), rxn1*k1*k2*rxn3, 1, 'lr');
%                 [vn,evn]=eigs(@(vec)(mv_t(mv(vec))), rxn1*k1*k2*rxn3, 1, 'sr');
%                 somedata{2}(i, swp) = sqrt(ev1/evn);
%                 fprintf('i=%d, cond(B): %g\n', i, somedata{2}(i, swp));

                if (strcmp(local_format, 'full'))
                    [sol_new,flg] = gmres(mv, rhs, nrestart, real_tol, 2, [], [], sol_prev);
                    %[dsol,flg]=gmres(mv, rhs-mv(sol_prev), nrestart, 1.0/8, 2, [], [], zeros(size(sol_prev)));
                    %sol_new=sol_prev+dsol;
                    res_new=norm(mv(sol_new)-rhs)/normf;
                    conv_factor=(res_new/res_prev);
                    if (res_new*(conv_factor)>real_tol && use_self_prec && strcmp(MatVec, 'bfun2')) % we need a prec.
                        if (strcmp(local_prec, 'selfprec'))
                            iB=tt_minres_selfprec(B, prec_tol, prec_compr, prec_iters, 'right');

                            resid = rhs-mv(sol_new);
                            [dsol,flg] = gmres(@(vec)bfun2(B, bfun2(iB,vec,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3),...
                                rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3), resid, nrestart, real_tol/res_new, gmres_iters);
                            dsol = bfun2(iB,dsol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
                            sol = sol_new+dsol;

                        end;
                        if (strcmp(local_prec, 'als'))
                            sol = als_solve_rx_2(B, rhs, real_tol, [], sol_new);
                        end;
                    else
                        [sol,flg] = gmres(mv, rhs, nrestart, real_tol, gmres_iters, [], [], sol_new);
                    end;
                    if (flg>0)
                        fprintf('-warn- gmres did not converge at block %d\n', i);
                    end;                    

                    res=mv(sol);
                    res_true = norm(res-rhs)/normf;
                else
                    sol_new = tt_gmres(mv, rhs, real_tol, 2, nrestart, real_tol, real_tol, [], [], [], sol_prev);
                    res_new=tt_dist3(mv(sol_new,[],[]),rhs)/normf;
                    conv_factor=(res_new/res_prev);
                    if (res_new*conv_factor>real_tol && use_self_prec && strcmp(MatVec, 'bfun2')) % we need a prec.
%                         if (strcmp(local_prec, 'selfprec'))
                            iB=tt_minres_selfprec(B, prec_tol, prec_compr, prec_iters, 'right');

%                             resid = tt_add(rhs, tt_scal(mv(sol_new,[],[]), -1));
%                             resid = tt_compr2(resid, real_tol);

                            sol = tt_gmres(@(vec,eps,mr)bfun3(B, vec, eps, mr), rhs, real_tol, gmres_iters, nrestart, real_tol, real_tol, @(vec,eps,mr)bfun3(iB, vec, eps, mr), [], [], sol_new);
%                             dsol = bfun3(iB,dsol,real_tol);
%                             sol = tt_add(sol_new,dsol);
%                             sol = tt_compr2(sol, real_tol);
%                         end;
                    else
                        sol = tt_gmres(mv, rhs, real_tol, gmres_iters, nrestart, real_tol, real_tol, [], [], [], sol_new);
                    end;
                    res=mv(sol,[],[]);
                    res_true = tt_dist3(res,rhs)/normf;
                end;
            end;

            if (strcmp(local_format, 'full'))
                dx(i) = norm(sol-sol_prev,'fro')/norm(sol_prev,'fro');
            else
                dx(i) = tt_dist3(sol, sol_prev)/exp(0.5*tt_dot2(sol,sol));
            end;
        else
            res_true = res_prev;
            dx(i)=0;
	    sol = sol_prev;
        end;

        if (verb>1)
            fprintf('=dmrg_solve2= Sweep %d, block %d, res_true = %3.3e\n', swp, i, res_true);
        end;
        if ((res_true>res_prev/trunc_to_true))&&(res_true>real_tol)&&(~last_sweep)
	    fprintf('--warn-- the residual damp by gmres was smaller than in the truncation\n');
%             keyboard;
            sol = sol_prev;
            res_true = res_prev;
        end;

        if (res_prev>max_res)
            max_res = res_prev;
        end;

        if (verb>1)
        fprintf('=dmrg_solve2= Sweep %d, block %d, dx=%3.3e, res_prev = %3.3e\n', swp, i, dx(i), res_prev);
        end;
        if (strcmp(local_format, 'full'))
            nrmsol = norm(sol, 'fro');
        else
            nrmsol = exp(0.5*tt_dot2(sol,sol));
        end
        if (nrmsol==0)
            dx(i)=0;
        end;

        if (swp==1)
            dx_old(i)=dx(i);
        end;

        % The new core does not converge - increase rank
        if (dx(i)/dx_old(i)>top_conv)&&(dx(i)>eps/(d^d_pow_check))
            drank(i)=drank(i)+ddrank;
            dpows(i)=dpows(i)+ddpow;
        end;
        % The new core converges well - try to decrease rank
        if (dx(i)/dx_old(i)<bot_conv)||(dx(i)<eps/(d^d_pow_check))
            drank(i)=max(drank(i)-ddrank, min_drank);
            dpows(i)=max(dpows(i)-ddpow, min_dpow);
        end;
        % perform simple compression for the last sweep
        if (last_sweep)
            dpows(i)=min(0.5, min_dpow);
        end;

        if (res_prev>bs_treshold*max_res_old)&&(strcmp(local_format, 'full'))
            if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
                [u,s,v]=svd(sol-reshape(sol_prev,[rxm1*m1,m2*rxm3]),'econ');
            else
                if (~last_sweep)
                    sol=reshape(sol,[rxm1*m1,m2*rxm3]);
                    [u,s,v]=svd(sol,'econ');
                else
                    [x2,rv]=qr(x2.', 0);
                    x1 = x1*(rv.');
                    [u,s,v]=svd(x1, 'econ');
                    v = x2*v;
                end;
            end;
            s = diag(s);
            flm=norm(s);
            %Truncation block. We have to make it smarter by binary search
            r0 = 1; rM = min(size(s,1),rmax); r = round((r0+rM)/2);
            while (rM-r0>2)
                er0=norm(s(r+1:numel(s)));
                if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
                    sol = sol_prev+reshape(u(:,1:r)*diag(s(1:r))*(v(:,1:r))',rxm1*m1*m2*rxm3, 1);
                else
                    sol = reshape(u(:,1:r)*diag(s(1:r))*(v(:,1:r))',rxm1*m1*m2*rxm3, 1);
                end;
                if (strcmp(MatVec,'full')||strcmp(MatVec,'half-full'))
                    resid = norm(B*sol-rhs)/norm(rhs);
                else
                    resid = norm(bfun2(B,sol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)-rhs)/norm(rhs);
                end;
                %             if ( verb>1 )
                %             fprintf('=dmrg_solve2= sweep %d, block %d, r=%d, resid=%g, er0=%g, MatVec=%s, rB=%d\n', swp, i, r, resid, er0/flm, MatVec, rB);
                %             end
                if ((resid<max(res_true*trunc_to_true, eps/(d^dpows(i)))) ) %Value of the rank is OK
                    rM = r-1;
                    r = round((r0+rM)/2);
                else %Is not OK.
                    r0 = r;
                    r = round((r0+rM)/2);
                end;
            end
            r = r0;
            % Line search - if the rank is underestimated
            cursol = cell(2,1);
            cursol{1}=u(:,1:r);
            cursol{2}=conj(v(:,1:r))*diag(s(1:r));
            if (strcmp(MatVec,'full')||strcmp(MatVec,'half-full'))
                resid = B*full(tt_tensor(cursol),rxm1*m1*m2*rxm3)-rhs;
                %resid = B*reshape(tt_to_full(cursol), rxm1*m1*m2*rxm3, 1)-rhs;                
            else
                resid = full(tt_tensor(tt_mv(B,cursol)),rxm1*m1*m2*rxm3)-rhs;
                %resid = reshape(tt_to_full(tt_mv(B,cursol)), rxm1*m1*m2*rxm3, 1)-rhs;
            end;
            while (r<min(size(s,1), rmax))
                r=r+1;
                er0=norm(s(r+1:numel(s)));
                cursol{1}=u(:,r);
                cursol{2}=conj(v(:,r))*s(r);
                if (strcmp(MatVec,'full')||strcmp(MatVec,'half-full'))
                    resid = B*full(tt_tensor(cursol),rxm1*m1*m2*rxm3)+resid;
                    %resid = B*reshape(tt_to_full(cursol), rxm1*m1*m2*rxm3, 1) + resid;

                else
                    resid = full(tt_tensor(tt_mv(B,cursol)),rxm1*m1*m2*rxm3)+resid;
                    %resid = reshape(tt_to_full(tt_mv(B,cursol)), rxm1*m1*m2*rxm3, 1)+resid;

                end;
                normres = norm(resid)/norm(rhs);
                %             if ( verb>1 )
                %                 fprintf('=dmrg_solve2= sweep %d, block %d, r=%d, resid=%g, er0=%g, MatVec=%s, rB=%d\n', swp, i, r, normres, er0/flm, MatVec, rB);
                %             end
                if ((normres<max(res_true*trunc_to_true, eps/(d^dpows(i)))) ) %Value of the rank is OK
                    break;
                end;
            end;

            if (~last_sweep)
                r = r+drank(i); % we want even larger ranks
            end;

            v = conj(v);
        else
            if (strcmp(local_format, 'tt'))
                x1 = sol{1};
                x2 = sol{2}.';
            end;
            % We do not have to decimate the whole supercore,
            % only one of factors, as we have the previous solution
            [v,rv]=qr(x2.',0); % size m2*rxm3, rxm2' - rxm2',rxm2
            r = size(v,2);
% %             rxm2 = size(x2,2);
            u = x1*rv.';
            s = ones(r,1);
%
%             [u,s,v]=svd(x1, 'econ');
%             v = x2*conj(v);
%             s = diag(s);
%             flm=norm(s);
%             %Truncation block. We have to make it smarter by binary search
%             r0 = 1; rM = min(size(s,1),rmax); r = round((r0+rM)/2);
%             while (rM-r0>2)
%                 er0=norm(s(r+1:numel(s)));
%                 if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
%                     sol = sol_prev+reshape(u(:,1:r)*diag(s(1:r))*(v(:,1:r)).',rxm1*m1*m2*rxm3, 1);
%                 else
%                     sol = reshape(u(:,1:r)*diag(s(1:r))*(v(:,1:r)).',rxm1*m1*m2*rxm3, 1);
%                 end;
%                 if (strcmp(MatVec,'full'))
%                     resid = norm(B*sol-rhs)/norm(rhs);
%                 else
%                     resid = norm(bfun2(B,sol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)-rhs)/norm(rhs);
%                 end;
%                 if ((resid<max(res_true*trunc_to_true, eps/(d^dpows(i)))) ) %Value of the rank is OK
%                     rM = r-1;
%                     r = round((r0+rM)/2);
%                 else %Is not OK.
%                     r0 = r;
%                     r = round((r0+rM)/2);
%                 end;
%             end
%             r = r0;
%             % Line search - if the rank is underestimated
%             cursol = cell(2,1);
%             cursol{1}=u(:,1:r);
%             cursol{2}=v(:,1:r)*diag(s(1:r));
%             if (strcmp(MatVec,'full'))
%                 resid = B*reshape(tt_to_full(cursol), rxm1*m1*m2*rxm3, 1)-rhs;
%             else
%                 resid = reshape(tt_to_full(tt_mv(B,cursol)), rxm1*m1*m2*rxm3, 1)-rhs;
%             end;
%             while (r<min(size(s,1), rmax))
%                 r=r+1;
%                 er0=norm(s(r+1:numel(s)));
%                 cursol{1}=u(:,r);
%                 cursol{2}=v(:,r)*s(r);
%                 if (strcmp(MatVec,'full'))
%                     resid = resid + B*reshape(tt_to_full(cursol), rxm1*m1*m2*rxm3, 1);
%                 else
%                     resid = resid + reshape(tt_to_full(tt_mv(B,cursol)), rxm1*m1*m2*rxm3, 1);
%                 end;
%                 normres = norm(resid)/norm(rhs);
%                 if ((normres<max(res_true*trunc_to_true, eps/(d^dpows(i)))) ) %Value of the rank is OK
%                     break;
%                 end;
%             end;
        end;
        r = min(r, max(size(s))); % but not too large
        r = min(r,rmax);

        v = v(:,1:r);
        u = u(:,1:r)*diag(s(1:r));

        if ( verb>1 )
            fprintf('=dmrg_solve2= sweep %d, block %d, r=%d, resid=%g, er0=%g, MatVec=%s, rB=%d\n', swp, i, r, normres, er0/flm, MatVec, rB);
        end

        % Keep rank increasing for several iterations
        % It helps for problems with hundred dimensions
%         if (mod(swp,dropsweeps)~=0)&&(dropflag==0)
%             r = max(r, ryold);
%         end;
%         if (verb>1)
%             fprintf('sweep %d, block %d, rank: %d, drop: %d\n', swp, i, r, dropflag);
%         end;
%         if (dropflag==1)&&(i==2)
%             dropflag=0;
%         end;

        % random kick %This is a production code, sir!
        %Replace by new stuff

        if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
            u = [x1, u];
            v = [x2.', v];
            [v,rv]=qr(v,0);
            u = u*(rv.');
            r = size(v,2);
        else
            if (~last_sweep)
                vr=randn(size(v,1),kickrank);
%                 v=reort(v,vr);
                [v,rr]=qr([v,vr], 0);
%                 radd=size(v,2)-r;
                radd=kickrank;
                if ( radd > 0 )
                    ur=zeros(size(u,1),radd);
                    u=[u,ur]*(rr.');
                end
                r = size(u,2);
%                 r=r+radd;
            end;
        end;

        %v = [v, randn(size(v,1), kickrank)];
        %u = [u, zeros(size(u,1), kickrank)];
        %[v,rv] = qr(v,0);
        %r = size(v,2);
        %u = u*(rv.');

        x{i}=permute(reshape(v, m2, rxm3, r), [1 3 2]);
        x{i-1}=permute(reshape(u, rxm1, m1, r), [2 1 3]);
        rxm2=r; rxn2=r;

        if (verb>2)
            % check the residual
            n1 = size(A{1},1);
            ra1 = size(A{1},4);
            rx1 = size(x{1},3);
            ry1 = size(y{1},3);
            A{1}=squeeze(A{1});
            x{1}=squeeze(x{1});
            y{1}=squeeze(y{1});
            true_resid = tt_mv(A,x);
            true_resid = tt_dist3(true_resid, y)/sqrt(tt_dot(y,y));
            fprintf('=dmrg_solve2= Sweep %d, block %d, true_resid: %3.3e\n', swp, i, true_resid);
            A{1}=reshape(A{1}, n1, n1, 1, ra1);
            x{1}=reshape(x{1}, n1, 1, rx1);
            y{1}=reshape(y{1}, n1, 1, ry1);
        end;



        phAold = reshape(phAold, rxn3*rp3*ra3, rxm3);
        x2 = reshape(x{i}, m2*rxm2, rxm3);
        phAold = phAold*(x2.'); % size rxn3*rp3*ra3*m2*rxm2
        phAold = reshape(phAold, rxn3, rp3, ra3, m2, rxm2);
        phAold = permute(phAold, [1 2 5 4 3]);
        phAold = reshape(phAold, rxn3*rp3*rxm2, m2*ra3);
        a2 = reshape(permute(A{i}, [2 4 1 3]), m2*ra3, n2*ra2);
        phAold=phAold*a2; % size rxn3*rp3*rxm2*n2*ra2
        phAold=reshape(phAold, rxn3, rp3, rxm2, n2, ra2);
        phAold = permute(phAold, [1 3 5 4 2]);
        phAold = reshape(phAold, rxn3*rxm2*ra2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        phAold=phAold*p2; % size rxn3*rxm2*ra2*k2*rp2
        phAold=reshape(phAold, rxn3,rxm2,ra2,k2,rp2);
        phAold = permute(phAold, [2 3 5 1 4]);
        phAold = reshape(phAold, rxm2*ra2*rp2, rxn3*k2);
        x2 = reshape(permute(x{i}, [3 1 2]), rxn3*k2, rxn2);
        phAold = phAold*conj(x2); % size rxm2*ra2*rp2*rxn2 <-- cplx conjugate!
        phAold = permute(reshape(phAold, rxm2, ra2, rp2, rxn2), [4 3 2 1]);

        phyold = reshape(phyold, rxn3*rp3, ry3);
        y2 = reshape(y{i}, n2*ry2, ry3);
        phyold = phyold*(y2.'); % size rxn3*rp3*n2*ry2
        phyold = reshape(phyold, rxn3, rp3, n2, ry2);
        phyold = permute(phyold, [1 4 3 2]);
        phyold = reshape(phyold, rxn3*ry2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        phyold = phyold*p2; % size rxn3*ry2*k2*rp2
        phyold = reshape(phyold, rxn3, ry2, k2, rp2);
        phyold = permute(phyold, [4 2 1 3]);
        phyold = reshape(phyold, rp2*ry2, rxn3*k2);
        x2 = reshape(permute(x{i}, [3 1 2]), rxn3*k2, rxn2);
        phyold = phyold*conj(x2); % size rp2*ry2*rxn2 <-- cplx conjugate!
        phyold = permute(reshape(phyold, rp2, ry2, rxn2), [3 1 2]);
    end;

%     sol_hist{3}=sol_hist{2};
%     sol_hist{2}=sol_hist{1};
%     sol_hist{1}=x;

    max_res_old = max_res;
    if (mod(swp,100)==0)
        max_res_old = 0;
    end;

%     if (max_res>max_res_old)
%         x = sol_hist{3};
%     else
%         max_res_old = max_res;
%     end;

%     if (max_res<tol*2*sqrt(d-1))
%         dropflag = 1;
%     end;
%     if (mod(swp,chksweeps)==0)||(swp==1)
%         x{1}=reshape(x{1}, size(x{1},1), size(x{1},3));
%         reschk = norm(tt_tensor(x)-tt_tensor(x_prev))/sqrt(tt_dot(x,x));
%         x_prev = x;
%         x{1}=reshape(x{1}, size(x{1},1),1, size(x{1},2));
%     end;
    if (verb>0)
        erank=0; sumn=0;
        for i=1:d
            erank = erank+size(x{i},1)*size(x{i},2)*size(x{i},3);
            sumn = sumn+size(x{i},1);
        end;
        erank = sqrt(erank/sumn);
        
%        A{1} = reshape(A{1}, size(A{1},1), size(A{1},2), size(A{1},4));
%        y{1} = reshape(y{1}, size(y{1},1), size(y{1},3));
%        x{1} = reshape(x{1}, size(x{1},1), size(x{1},3));
%        res_real = tt_dist3(tt_mv(A, x), y)/sqrt(tt_dot(y,y));
%        A{1} = reshape(A{1}, size(A{1},1), size(A{1},2), 1, size(A{1},3));
%        y{1} = reshape(y{1}, size(y{1},1), 1, size(y{1},2));
%        x{1} = reshape(x{1}, size(x{1},1), 1, size(x{1},2));        
%        somedata{4}(swp) = res_real;
        
%         fprintf('===Sweep %d, res_%d: %3.3e, drop_next: %d, dx_max: %3.3e, res_max: %3.3e\n', swp, chksweeps,0, dropflag, dx_max, max_res);
        fprintf('=dmrg_solve2= Sweep %d, dx_max: %3.3e, res_max: %3.3e, erank: %g\n', swp, max(dx), max_res, erank);
    end;
    if (last_sweep)
        break;
    end;
    if (max_res<tol/(d^d_pow_check))||(swp==nswp-1)
        last_sweep=true;
%         break;
    end;

    dx_old = dx;
%     if (verb>0)
%         fprintf('-=-=-=-=-= dx_max = %3.3e, res_max = %3.3e\n', dx_max, max_res);
%     end;
%     if (dx_max<tol*2)
%     if (max_res<tol*2)
%         break;
%     end;
%     keyboard;
end;

x{1}=reshape(x{1}, size(x{1},1), size(x{1},3));

% x = tt_compr2(x, eps, rmax);

if (input_is_tt_tensor)
  x=tt_tensor(x);
end

if (nargout>1)
    somedata{1} = swp;
end;

end

function [y]=bfun2(B, x, rxm1, m1, m2, rxm3, rxn1, k1, k2, rxn3)
% Computes (B{1} \otimes B{2})x
% B{1} is of sizes rxn1*k1, rxm1*m1, rB
% B{2} is of sizes k2*rxn3, m2*rxm3, rB
rB=size(B{1},3);
x = reshape(x, rxm1*m1, m2*rxm3);
B1 = permute(B{1}, [3 1 2]);
B1 = reshape(B1, rB*rxn1*k1, rxm1*m1);
y = B1*x; % size rB*rxn1*k1,m2*rxm3
y = reshape(y, rB, rxn1*k1, m2*rxm3);
y = permute(y, [3 1 2]);
y = reshape(y, m2*rxm3*rB, rxn1*k1);
B2 = reshape(B{2}, k2*rxn3, m2*rxm3*rB);
y = B2*y; % size k2*rxn3,rxn1*k1
y = reshape(y.', rxn1*k1*k2*rxn3, 1);
end


function [y]=bfun2_t(B, x, rxm1, m1, m2, rxm3, rxn1, k1, k2, rxn3)
% Computes (B{1}' \otimes B{2}')x
% B{1} is of sizes rxn1*k1, rxm1*m1, rB
% B{2} is of sizes k2*rxn3, m2*rxm3, rB
rB=size(B{1},3);
x = reshape(x, rxm1*m1, m2*rxm3);
B1 = permute(B{1}, [3 2 1]);  % HERE
B1 = reshape(B1, rB*rxn1*k1, rxm1*m1);
y = B1*x; % size rB*rxn1*k1,m2*rxm3
y = reshape(y, rB, rxn1*k1, m2*rxm3);
y = permute(y, [3 1 2]);
y = reshape(y, m2*rxm3*rB, rxn1*k1);
B2 = reshape(permute(B{2}, [2,1,3]), k2*rxn3, m2*rxm3*rB);
y = B2*y; % size k2*rxn3,rxn1*k1
y = reshape(y.', rxn1*k1*k2*rxn3, 1);
end


function [y] = bfun3(A, x, eps, mr)
% For the 2d--TT MatVec
y = tt_mv(A, x);
if (nargin<4)
    mr = [];
end;
if (nargin>2)&&(~isempty(eps))
    y = tt_compr2(y, eps, mr);
end;

end

