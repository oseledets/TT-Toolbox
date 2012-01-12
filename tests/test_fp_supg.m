dconf = 1;
dpx = 2;
d0x = 8;

a = 10;
beta = 0.5;

eps = 1e-8;
tol = 1e-6;

timescheme = 'ts';
Nt = 1024;
T = 10;

%  timescheme = 'block';

% diffusion tensor - in conf vars
Dc = spdiags(ones(dconf,1)*[-1,2,-1], [-1,0,1], dconf,dconf);
Dc = full(Dc);
D = kron(Dc, eye(dpx)); % for all variables

v1 = @(x)(x); % elementary velocity
vp1 = @(x)(1+0*x); % elementary velocity gradient

% Ext. flow - sparse format
Ki = [1];
Kj = [2];
Kel = beta;


% Grid
N = 2^d0x;
h = 2*a/(N+1);
x = (-a+h:h:a-h)'; % comp. domain
xb = (-a:h:a)';  % + boundaries

% Elementary matrices
M = full(spdiags(ones(N,1)*[1/6, 4/6, 1/6]*h, [-1,0,1], N, N));
S1 = full(spdiags(ones(N,1)*[-0.5, 1, 0.5], [-1,0,1], N, N));
S1const = full(spdiags(ones(N,1)*[-0.5, 0, 0.5], [-1,0,1], N, N));
S2 = full(spdiags(ones(N,1)*[-1, 2, -1]/h, [-1,0,1], N, N));

% Full-boundary and inner coordinates
Xb = cell(dpx*dconf, 1);
X = cell(dpx*dconf, 1);
for i=1:dpx*dconf
    for j=1:dpx*dconf
	if (i==j)
	    Xb{i} = kron(Xb{i}, tt_tensor(xb));
	    X{i} = kron(X{i}, tt_tensor(x));
	else
	    Xb{i} = kron(Xb{i}, tt_tensor(ones(size(xb))));
	    X{i} = kron(X{i}, tt_tensor(ones(size(x))));
	end;
    end;
end;

% Prepare velocities
% First, elementary
Vb_el = cell(dpx*dconf, 1);
V_el = cell(dpx*dconf, 1);
for i=1:dpx*dconf
    for j=1:dpx*dconf
	if (i==j)
	    curVb = v1(xb);
	    curVb(isnan(curVb))=0;
	    curVb(isinf(curVb))=0;
	    curVb = tt_tensor(curVb);
	    curV = v1(x);
	    curV(isnan(curV))=0;
	    curV(isinf(curV))=0;
	    curV = tt_tensor(curV);	    
	else
	    curVb = tt_tensor(ones(size(xb)));
	    curV = tt_tensor(ones(size(x)));
	end;
	Vb_el{i} = kron(Vb_el{i}, curVb);
	V_el{i} = kron(V_el{i}, curV);
    end;
end;
% Now, real velocities to stuff into the equation
Vb = cell(dpx*dconf, 1);
V = cell(dpx*dconf, 1);
for i=1:dconf
    for j=1:dpx
	% Rouse-components
	for k=1:dconf
	    Vb{j+(i-1)*dpx} = Vb{j+(i-1)*dpx}-0.25*Dc(i,k)*Vb_el{j+(k-1)*dpx};
	    V{j+(i-1)*dpx} = V{j+(i-1)*dpx}-0.25*Dc(i,k)*V_el{j+(k-1)*dpx};
	end;
	% Ext. flow comps.
	extind = find(Ki==j);
	if (~isempty(extind))
	    Vb{j+(i-1)*dpx} = Vb{j+(i-1)*dpx} + Kel(extind)*Xb{Kj(extind)+(i-1)*dpx};
	    V{j+(i-1)*dpx} = V{j+(i-1)*dpx} + Kel(extind)*X{Kj(extind)+(i-1)*dpx};
	end;
	Vb{j+(i-1)*dpx} = round(Vb{j+(i-1)*dpx}, eps);
	V{j+(i-1)*dpx} = round(V{j+(i-1)*dpx}, eps);
    end;
end;

% Velocity gradients
% elementary
Vpb_el = cell(dpx*dconf,1);
for i=1:dpx*dconf
    for j=1:dpx*dconf
	if (i==j)
	    curV = vp1(xb);
	    curV(isnan(curV))=0;
	    curV(isinf(curV))=0;
	    curV = tt_tensor(curV);	    
	else
	    curV = tt_tensor(ones(size(xb)));
	end;
	Vpb_el{i} = kron(Vpb_el{i}, curV);
    end;
end;
% Real
Vpb = cell(dpx*dconf, 1);
for i=1:dconf
    for j=1:dpx
	% Rouse-components
	for k=1:dconf
	    Vpb{j+(i-1)*dpx} = Vpb{j+(i-1)*dpx}-0.25*Dc(i,k)*Vpb_el{j+(k-1)*dpx};
	end;
	% Ext. flow comps.
	extind = find(Ki==j);
	if (~isempty(extind))
	    % Gradient of Kx here
	    Vpb{j+(i-1)*dpx} = Vpb{j+(i-1)*dpx} + Kel(extind)*tt_tensor(tt_ones(dpx*dconf,N+2));
	end;
	Vpb{j+(i-1)*dpx} = round(Vpb{j+(i-1)*dpx}, eps);
    end;
end;
% The whole divergence - our reaction coefficient
Vp = [];
for i=1:dpx*dconf
    Vp = Vp+Vpb{i};
end;
Vp = round(Vp, eps);

%  keyboard;
% Now, generate the coefs in half-int points
for i=1:dpx*dconf
    Vi = Vb{i};
    for j=1:dpx*dconf
	curV = Vi{j};
	r1 = size(curV,1); r2 = size(curV,3);
	curV = permute(curV, [2, 1, 3]);
	curV = reshape(curV, N+2, r1*r2);
	curV = (curV(1:N+1, :)+curV(2:N+2, :))*0.5;
	curV = reshape(curV, N+1, r1, r2);
	Vi{j} = permute(curV, [2, 1, 3]);
    end;
    Vb{i} = Vi;
end;
for j=1:dpx*dconf
    curV = Vp{j};
    r1 = size(curV,1); r2 = size(curV,3);
    curV = permute(curV, [2, 1, 3]);
    curV = reshape(curV, N+2, r1*r2);
    curV = (curV(1:N+1, :)+curV(2:N+2, :))*0.5;
    curV = reshape(curV, N+1, r1, r2);
    Vp{j} = permute(curV, [2, 1, 3]);    
end;

%  keyboard;
% Finally, we are ready to generate some matrices
% Second Partial derivatives
Stiff2 = [];
for i=1:dpx*dconf
    for j=1:dpx*dconf
	if (D(i,j)~=0)
	    curS = [];
	    for k=1:dpx*dconf
		if (k==i)&&(k==j)
		    curS = kron(curS, tt_matrix(S2));
		elseif (k==i)&&(k~=j)
		    curS = kron(curS, tt_matrix(S1const'));
		elseif (k==j)&&(k~=i)
		    curS = kron(curS, tt_matrix(S1const));
		else
		    curS = kron(curS, tt_matrix(M));
		end;
	    end;
	    Stiff2 = Stiff2 + 0.25*D(i,j)*curS;
	    Stiff2 = round(Stiff2, eps);
	end;
    end;
end;

% Advection
Stiff1 = [];
for i=1:dpx*dconf
    Vi = Vb{i};
    curS = [];
    for j=1:dpx*dconf
	curV = Vi{j};
	r1 = size(curV,1); r2 = size(curV,3);
	curV = permute(curV, [2, 1, 3]);
	curV = reshape(curV, N+1, r1*r2);
	Sbl = zeros(N,N,r1*r2);
	for k=1:r1*r2
	    vdiag = zeros(N,3);
	    vdiag(1:N-1,1) = curV(2:N,k);
	    vdiag(2:N,3) = curV(2:N,k);	    
	    if (i==j)
		vdiag(1:N,2) = (curV(1:N,k)-curV(2:N+1,k))*0.5;
		vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
		vdiag = vdiag.*S1;
	    else
		vdiag(1:N,2) = (curV(1:N,k)+curV(2:N+1,k))*0.5;
		vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));    
		vdiag = vdiag.*M;
	    end;
	    Sbl(:,:,k) = vdiag;
	end;
	Sbl = reshape(Sbl, N, N, r1, r2);
	Sbl = permute(Sbl, [3, 1, 2, 4]);
%  	Sbl = reshape(Sbl, r1*N*N*r2, 1);
	Sblm = tt_matrix;
	Sblm.d = 1; Sblm.ps=ones(2,1); Sblm.r = ones(2,1);
	Sblm{1} = Sbl;
	Sbl = Sblm;
%  	Sbl = tt_matrix(Sbl);
%  	Sbl.n = N; Sbl.m = N;
%  	Sbl.r = [r1; r2];
	curS = kron(curS, Sbl);
    end;
    Stiff1 = Stiff1 + curS;
    Stiff1 = round(Stiff1, eps);
end;
%  keyboard;

% reaction
Stiff0 = [];
for j=1:dpx*dconf
    curV = Vp{j};
    r1 = size(curV,1); r2 = size(curV,3);
    curV = permute(curV, [2, 1, 3]);
    curV = reshape(curV, N+1, r1*r2);
    Sbl = zeros(N,N,r1*r2);
    for k=1:r1*r2
	sdiag = zeros(N,3);
	sdiag(1:N-1,1) = curV(2:N,k);
	sdiag(1:N,2) = (curV(1:N,k)+curV(2:N+1,k))*0.5;
	sdiag(2:N,3) = curV(2:N,k);
	sdiag = full(spdiags(sdiag, [-1, 0, 1], N, N));    
	sdiag = sdiag.*M;
	Sbl(:,:,k) = sdiag;
    end;
    Sbl = reshape(Sbl, N, N, r1, r2);
    Sbl = permute(Sbl, [3, 1, 2, 4]);
    Sblm = tt_matrix;
    Sblm.d = 1; Sblm.ps=ones(2,1); Sblm.r = ones(2,1);
    Sblm{1} = Sbl;
    Sbl = Sblm;    
%      Sbl = reshape(Sbl, r1*N, N*r2);
%      Sbl = tt_matrix(Sbl);
%      Sbl.n = N; Sbl.m = N;
%      Sbl.r = [r1; r2];
    Stiff0 = kron(Stiff0, Sbl);
end;

% The basic stiffness matrix is ready
Ax = Stiff2 + Stiff1 + Stiff0;
Ax = round(Ax, eps);

% The basic mass matrix - for time derivative
Mass = [];
for i=1:dpx*dconf
    Mass = kron(Mass, tt_matrix(M));
end;

% Now, the most serios
%  SUPG part
% tau
tau = cell(dpx*dconf, 1);
for i=1:dpx*dconf
    Vi = V{i};
    Vi = tt_reshape(Vi, 2*ones(d0x*dpx*dconf, 1), eps);
    tau{i} = funcrs2(Vi, @(v)supg_tau(v, 0.25*D(i,i), h), eps, Vi, 20);
%      tau{i} = funcrs2(Vi, @(v)(h*0.5*(coth(max(abs(v*h),1e-8))-1./max(abs(v*h),1e-8))./max(abs(v),1e-8)), eps, Vi, 20);
    tau{i} = tt_reshape(tau{i}, N*ones(dpx*dconf, 1));
    tau{i} = diag(tau{i});
end;

% Convection-Convection
SUPG_VV = [];
for i=1:dpx*dconf
    for j=1:dpx*dconf
	Vij = Vb{i}.*Vb{j};
	Vij = round(Vij, eps);
	curS = [];
	for k=1:dpx*dconf
	    curV = Vij{k};
	    r1 = size(curV,1); r2 = size(curV,3);
	    curV = permute(curV, [2, 1, 3]);
	    curV = reshape(curV, N+1, r1*r2);
	    Sbl = zeros(N,N,r1*r2);
	    for r=1:r1*r2
		vdiag = zeros(N,3);
		vdiag(1:N-1,1) = curV(2:N,r);
		vdiag(2:N,3) = curV(2:N,r);
		if (k==i)&&(k==j)
		    vdiag(1:N,2) = (curV(1:N,r)+curV(2:N+1,r))*0.5;
		    vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
		    vdiag = vdiag.*S2;
		elseif (k==i)&&(k~=j)
		    vdiag(1:N,2) = (curV(1:N,r)-curV(2:N+1,r))*0.5;
		    vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));		
		    vdiag = vdiag.*(S1');
		elseif (k==j)&&(k~=i)
		    vdiag(1:N,2) = (curV(1:N,r)-curV(2:N+1,r))*0.5;
		    vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));		
		    vdiag = vdiag.*S1;
		else
		    vdiag(1:N,2) = (curV(1:N,r)+curV(2:N+1,r))*0.5;
		    vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
		    vdiag = vdiag.*M;
		end;
		Sbl(:,:,r) = vdiag;
	    end;
	    Sbl = reshape(Sbl, N, N, r1, r2);
	    Sbl = permute(Sbl, [3, 1, 2, 4]);
	    Sblm = tt_matrix;
	    Sblm.d = 1; Sblm.ps=ones(2,1); Sblm.r = ones(2,1);
	    Sblm{1} = Sbl;
	    Sbl = Sblm;
%  	    Sbl = reshape(Sbl, r1*N, N*r2);
%  	    Sbl = tt_matrix(Sbl);
%  	    Sbl.n = N; Sbl.m = N;
%  	    Sbl.r = [r1; r2];
	    curS = kron(curS, Sbl);
	end;
	SUPG_VV = SUPG_VV + tau{i}*curS;
	SUPG_VV = round(SUPG_VV, eps);
    end;
end;

% Convection-reaction
SUPG_VS = [];
for i=1:dpx*dconf
    VSi = Vb{i};
    VSi = VSi.*Vp;
    VSi = round(VSi, eps);
    curS = [];
    for j=1:dpx*dconf
	curV = VSi{j};
	r1 = size(curV,1); r2 = size(curV,3);
	curV = permute(curV, [2, 1, 3]);
	curV = reshape(curV, N+1, r1*r2);
	Sbl = zeros(N,N,r1*r2);
	for k=1:r1*r2
	    vdiag = zeros(N,3);
	    vdiag(1:N-1,1) = curV(2:N,k);
	    vdiag(2:N,3) = curV(2:N,k);
	    if (i==j)
		vdiag(1:N,2) = (curV(1:N,k)-curV(2:N+1,k))*0.5;
		vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
		vdiag = vdiag.*(S1');
	    else
		vdiag(1:N,2) = (curV(1:N,k)+curV(2:N+1,k))*0.5;
		vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
		vdiag = vdiag.*M;
	    end;
	    Sbl(:,:,k) = vdiag;
	end;
	Sbl = reshape(Sbl, N, N, r1, r2);
	Sbl = permute(Sbl, [3, 1, 2, 4]);
	Sblm = tt_matrix;
	Sblm.d = 1; Sblm.ps=ones(2,1); Sblm.r = ones(2,1);
	Sblm{1} = Sbl;
	Sbl = Sblm;
%  	Sbl = reshape(Sbl, r1*N, N*r2);
%  	Sbl = tt_matrix(Sbl);
%  	Sbl.n = N; Sbl.m = N;
%  	Sbl.r = [r1; r2];
	curS = kron(curS, Sbl);
    end;
    SUPG_VS = SUPG_VS + tau{i}*curS;
    SUPG_VS = round(SUPG_VS, eps);
end;

% Convection-time part
SUPG_VT = [];
for i=1:dpx*dconf
    VSi = Vb{i};
    curS = [];
    for j=1:dpx*dconf
	curV = VSi{j};
	r1 = size(curV,1); r2 = size(curV,3);
	curV = permute(curV, [2, 1, 3]);
	curV = reshape(curV, N+1, r1*r2);
	Sbl = zeros(N,N,r1*r2);
	for k=1:r1*r2
	    vdiag = zeros(N,3);
	    vdiag(1:N-1,1) = curV(2:N,k);
	    vdiag(2:N,3) = curV(2:N,k);
	    if (i==j)
		vdiag(1:N,2) = (curV(1:N,k)-curV(2:N+1,k))*0.5;
		vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
		vdiag = vdiag.*(S1');
	    else
		vdiag(1:N,2) = (curV(1:N,k)+curV(2:N+1,k))*0.5;
		vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
		vdiag = vdiag.*M;
	    end;
	    Sbl(:,:,k) = vdiag;
	end;
	Sbl = reshape(Sbl, N, N, r1, r2);
	Sbl = permute(Sbl, [3, 1, 2, 4]);
	Sblm = tt_matrix;
	Sblm.d = 1; Sblm.ps=ones(2,1); Sblm.r = ones(2,1);
	Sblm{1} = Sbl;
	Sbl = Sblm;
%  	Sbl = reshape(Sbl, r1*N, N*r2);
%  	Sbl = tt_matrix(Sbl);
%  	Sbl.n = N; Sbl.m = N;
%  	Sbl.r = [r1; r2];
	curS = kron(curS, Sbl);
    end;
    SUPG_VT = SUPG_VT + tau{i}*curS;
    SUPG_VT = round(SUPG_VT, eps);
end;

% Convection-Diffusion
SUPG_VD = [];
for i=1:dpx*dconf
    for j=1:dpx*dconf
	for k=1:dpx*dconf
	if ((i~=j)||(i~=k))&&(D(i,j)~=0)
	    Vk = Vb{k};
	    curS = [];
	    for m=1:dpx*dconf
		curV = Vk{m};
		r1 = size(curV,1); r2 = size(curV,3);
		curV = permute(curV, [2, 1, 3]);
		curV = reshape(curV, N+1, r1*r2);
		Sbl = zeros(N,N,r1*r2);
		for r=1:r1*r2
		    vdiag = zeros(N,3);
		    vdiag(1:N-1,1) = curV(2:N,r);
		    vdiag(2:N,3) = curV(2:N,r);
		    if ((m==i)&&(m==j)&&(m~=k))||((m~=i)&&(m==j)&&(m==k))||((m==i)&&(m~=j)&&(m==k))
			% S2 on any
			vdiag(1:N,2) = (curV(1:N,r)+curV(2:N+1,r))*0.5;
			vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
			vdiag = vdiag.*S2;
		    elseif ((m~=i)&&(m==j)&&(m~=k))
			% S1 on j - no transpose
			vdiag(1:N,2) = (curV(1:N,r)-curV(2:N+1,r))*0.5;
			vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));		
			vdiag = vdiag.*S1;
		    elseif ((m==i)&&(m~=j)&&(m~=k))||((m~=i)&&(m~=j)&&(m==k))
			% S1 on i or k - transpose
			vdiag(1:N,2) = (curV(1:N,r)-curV(2:N+1,r))*0.5;
			vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
			vdiag = vdiag.*(S1');
		    else
			% m does not touch current dimensions - M
			vdiag(1:N,2) = (curV(1:N,r)+curV(2:N+1,r))*0.5;
			vdiag = full(spdiags(vdiag, [-1, 0, 1], N, N));
			vdiag = vdiag.*M;
		    end;
		    Sbl(:,:,r) = vdiag;
		end;
		Sbl = reshape(Sbl, N, N, r1, r2);
		Sbl = permute(Sbl, [3, 1, 2, 4]);
		Sblm = tt_matrix;
		Sblm.d = 1; Sblm.ps=ones(2,1); Sblm.r = ones(2,1);
		Sblm{1} = Sbl;
		Sbl = Sblm;		
%  		Sbl = reshape(Sbl, r1*N, N*r2);
%  		Sbl = tt_matrix(Sbl);
%  		Sbl.n = N; Sbl.m = N;
%  		Sbl.r = [r1; r2];
		curS = kron(curS, Sbl);
	    end;
	    SUPG_VD = SUPG_VD + 0.25*D(i,j)*tau{k}*curS;
	    SUPG_VD = round(SUPG_VD, eps);
	end;
	end;
    end;
end;

% Spatial SUPG matrix
AxSUPG = SUPG_VV + SUPG_VS;
if (~isempty(SUPG_VD))
    AxSUPG = AxSUPG + SUPG_VD;
end;
AxSUPG = round(AxSUPG, eps);


% Convert all matrices to the Qtt_tucker
Ax = tt_to_qtt3(core(Ax), 1, eps);
Ax = tt_matrix(Ax);
Ax = matrix(qtt_tucker(Ax.tt, d0x*ones(dpx*dconf, 1), eps));
Mass = tt_to_qtt3(core(Mass), 1, eps);
Mass = tt_matrix(Mass);
Mass = matrix(qtt_tucker(Mass.tt, d0x*ones(dpx*dconf, 1), eps));
AxSUPG = tt_to_qtt3(core(AxSUPG), 1, eps);
AxSUPG = tt_matrix(AxSUPG);
AxSUPG = matrix(qtt_tucker(AxSUPG.tt, d0x*ones(dpx*dconf, 1), eps));
SUPG_VT = tt_to_qtt3(core(SUPG_VT), 1, eps);
SUPG_VT = tt_matrix(SUPG_VT);
SUPG_VT = matrix(qtt_tucker(SUPG_VT.tt, d0x*ones(dpx*dconf, 1), eps));


% Now, prepare initial guess, time matrices, etc
uSN = exp(-0.5*(x.^2));
uSN = tt_tensor(uSN);
uSN = tt_reshape(uSN, 2*ones(d0x, 1), eps);
u0 = [];
for i=1:dpx*dconf
    u0 = kron(u0, uSN);
end;
u0 = qtt_tucker(u0, d0x*ones(dpx*dconf, 1), eps);


if (strcmp(timescheme, 'ts'))
    tau = T/Nt;
    
    CNm = Mass/tau - Ax*0.5   + SUPG_VT/tau - AxSUPG*0.5;
    CNm = round(CNm, eps);
    CNp = Mass/tau + Ax*0.5   + SUPG_VT/tau + AxSUPG*0.5;
    CNp = round(CNp, eps);  
    
    u = u0;
    sumtime = 0;
    for t=1:Nt
	tic;
	rhs = mvrk(CNm, u, tol);
	u = dmrg_rake_solve2(CNp, rhs, tol, 'x0', u);
	curtime = toc;
	sumtime = sumtime+curtime;
	
	Au = norm(Ax*u)/norm(Mass*u);
	
	u2 = qtttucker_to_tt(u.tuck, u.core);
	ind = num2cell((2^(d0x-6):2^(d0x-6):2^d0x)'*ones(1,2), 1);
	contour(full(u2(ind), 2^6*[1,1]));
%  	plot(full(u, N));
	str = sprintf('time: %d [%g], rank: %d', t, t*tau, max(max(rank(u))));
	title(str);
	print('-depsc', 'imgout');
%  	mybatchprint;
	
	fprintf('ts: %d [%g], Au: %3.3e, curtime: %g, tot. time: %g\n', t, t*tau, Au, curtime, sumtime);
	pause(0.1);
    end;
end;
