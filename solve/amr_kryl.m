function [V,y]=amr_kryl(A, funA, x, tol, varargin)
%Approximate the matrix function by vector via the AMEn+Krylov iteration
%   [V,y]=amr_kryl(A, funA, x, tol, varargin)
%   Attempts to approximate the Krylov basis for the problem w = f(A)*x
%   with accuracy TOL using the AMEn+ALS iteration.
%   funA is a function of type @(A,x), and should return f(A)*x
%   Matrix A has to be given in the TT-format, right-hand side x should be
%   given in the TT-format also. 
%   Returns the Krylov subspace V in the TT format with the last rank
%   enumerating the Krylov vectors, and the Krylov coefficients y.
%   The actual vector f(A)*x is constructed as V*y.
%
%   Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following)
%   The list of option names and default values are:
%       o v0 - initial approximation to Krylov basis [x]
%       o nswp - maximal number of sweeps [20]
%       o m - maximal Krylov dimension [20]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o kickrank - compression rank of the next Krylov vector, 
%         i.e. enrichment size [3]
%
%       Example (matrix exponential):
%           d=6; f=2; tol = 1e-5;
%           mat=-tt_qlaplace_dd(d*ones(1,f)); %Laplace in the QTT-format
%           rhs=tt_ones(2,d*f); Right-hand side of all ones
%           [V,y] = amr_kryl(mat, @(A,x)(expm(A)*x), rhs, tol); %Krylov-exp
%           sol = V*y; % evaluate the solution
%
%********
%   For description of adaptive ALS please see
%   S. Dolgov, D. Savostyanov,
%   Alternating minimal energy methods for linear systems
%   in higher dimensions. Part I: SPD systems,
%   http://arxiv.org/abs/1301.6068
%
%   Part II with more relevant details will come soon...
%   Use {sergey.v.dolgov, drraug}@gmail.com for feedback
%********
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



nswp = 20;
kickrank = 3;
Mmax = 20;
verb = 1;

V = x;

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'm'
            Mmax=varargin{i+1};
        case 'v0'
            V=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
            
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

d = x.d;
n = x.n;
rx = x.r;
rv = V.r;

x = core2cell(x);
V = core2cell(V);

A = core2cell(A);

if (kickrank>0)
    z = tt_rand(n, d, kickrank);
    rz = z.r;
    z = core2cell(z);
    phizax = cell(d+1,1); phizax{1}=1; phizax{d+1}=1;
end;

phia = cell(d+1,1); phia{1}=1; phia{d+1}=1;
phix = cell(d+1,1); phix{1}=1; phix{d+1}=1;

% Initial ort
for i=1:d-1
    cr = reshape(V{i}, rv(i)*n(i), rv(i+1));
    [cr,R]=qr(cr, 0);
    cr2 = reshape(V{i+1}, rv(i+1), n(i+1)*rv(i+2));
    cr2 = R*cr2;
    rv(i+1) = size(cr, 2);
    V{i} = reshape(cr, rv(i), n(i), rv(i+1));
    V{i+1} = reshape(cr2, rv(i+1), n(i+1), rv(i+2));
    phia{i+1} = compute_next_Phi(phia{i}, V{i}, A{i}, V{i}, 'lr');
    phix{i+1} = compute_next_Phi(phix{i}, V{i}, [], x{i}, 'lr');
    
    if (kickrank>0)
        cr = reshape(z{i}, rz(i)*n(i), rz(i+1));
        [cr,R]=qr(cr, 0);
        cr2 = reshape(z{i+1}, rz(i+1), n(i+1)*rz(i+2));
        cr2 = R*cr2;
        rz(i+1) = size(cr, 2);
        z{i} = reshape(cr, rz(i), n(i), rz(i+1));
        z{i+1} = reshape(cr2, rz(i+1), n(i+1), rz(i+2));
        phizax{i+1} = compute_next_Phi(phizax{i}, z{i}, A{i}, V{i}, 'lr');
    end;
end;

i = d;
dir = -1;
swp = 1;
max_dx = 0;
b = rv(d+1); % Initially, the enum rank is in V{d}
V{d} = reshape(V{d}, rv(d)*n(d), b);
% Add space for a new Krylov vector
if ((b+1)<=Mmax)
    V{d} = [V{d}, zeros(rv(d)*n(d), 1)];
    b=b+1;
end;
V{d} = reshape(V{d}, rv(d), n(d), 1, b);
rv(d+1) = 1;
b_new = 1;

while (swp<=nswp)
    % Project the generating vector
    crx = reshape(x{i}, rx(i), n(i)*rx(i+1));
    crx = phix{i}*crx;
    crx = reshape(crx, rv(i)*n(i), rx(i+1));
    crx = crx*phix{i+1};
    crx = reshape(crx, rv(i)*n(i)*rv(i+1), 1);
    
    % Build projected Krylov basis
    crv = zeros(rv(i)*n(i)*rv(i+1), b);
    crv(:,1) = crx/norm(crx);
    for j=2:b
        w = bfun3(phia{i}, A{i}, phia{i+1}, crv(:,j-1));
        crv(:,j)=w/norm(w);
    end;
    [crv,R]=qr(crv, 0);
    % Caution: b could be reduced
    b = size(crv, 2);
    
    % Measure the weights using the actual matrix function
    w = bfun3(phia{i}, A{i}, phia{i+1}, crv);
    B = crv'*w;
    y = funA(B, eye(b,1)*(crv(:,1)'*crx));
    inew = find(abs(y)/abs(y(1))<tol, 1);
    if (isempty(inew))
        inew = b;
    end;
    b_new = max(b_new, inew);
    
    % Caution #2: in the previous iteration there could be different b
    V_prev = reshape(V{i}, rv(i)*n(i)*rv(i+1), size(V{i}, 4));
    [V_prev,R]=qr(V_prev,0);
    w = bfun3(phia{i}, A{i}, phia{i+1}, V_prev);
    y_prev = funA((V_prev'*w), (V_prev'*crx));
    dx = norm(crv*y-V_prev*y_prev)/norm(y);
    
    max_dx = max(max_dx, dx);
    
    % Truncation and enrichment
    if ((dir>0)&&(i<d))||((dir<0)&&(i>1))
        % Apply weights
        crv = crv.*(ones(rv(i)*n(i)*rv(i+1), 1)*abs(y.'));
        if (dir>0)
            crv = reshape(crv, rv(i)*n(i), rv(i+1)*b);
            [u,s,v]=svd(crv, 'econ');
            s = diag(s);
            r = my_chop2(s, tol*norm(s)/sqrt(d));
            u = u(:,1:r);
            v = conj(v(:,1:r))*diag(s(1:r));
            
            % Prepare enrichment, if needed
            if (kickrank>0)
                crv = u*v.';
                crv = reshape(crv, rv(i)*n(i)*rv(i+1), b);
                % Take next Krylov vector A*V_b into z
                crv = crv(:,b);
                % For updating z
                crz = bfun3(phizax{i}, A{i}, phizax{i+1}, crv);
                % For adding into solution
                crs = bfun3(phia{i}, A{i}, phizax{i+1}, crv);
                crs = reshape(crs, rv(i)*n(i), rz(i+1));
                u = [u,crs];
                [u,R]=qr(u, 0);
                v = [v, zeros(rv(i+1)*b, rz(i+1))];
                v = v*R.';
                r = size(u, 2);
            end;
            V{i} = reshape(u, rv(i), n(i), r);
            
            cr2 = reshape(V{i+1}, rv(i+1), n(i+1)*rv(i+2));
            v = reshape(v, rv(i+1), b*r);
            cr2 = v.'*cr2;
            cr2 = reshape(cr2, b, r, n(i+1), rv(i+2));
            V{i+1} = permute(cr2, [2,3,4,1]);
            
            rv(i+1) = r;
            
            phia{i+1} = compute_next_Phi(phia{i}, V{i}, A{i}, V{i}, 'lr');
            phix{i+1} = compute_next_Phi(phix{i}, V{i}, [], x{i}, 'lr');
            
            if (kickrank>0)
                crz = reshape(crz, rz(i)*n(i), rz(i+1));
                [crz,R]=qr(crz, 0);
                cr2 = reshape(z{i+1}, rz(i+1), n(i+1)*rz(i+2));
                cr2 = R*cr2;
                rz(i+1) = size(crz, 2);
                z{i} = reshape(crz, rz(i), n(i), rz(i+1));
                z{i+1} = reshape(cr2, rz(i+1), n(i+1), rz(i+2));
                
                phizax{i+1} = compute_next_Phi(phizax{i}, z{i}, A{i}, V{i}, 'lr');
            end;
        else       
            crv = reshape(crv.', b*rv(i), n(i)*rv(i+1));            
            [u,s,v]=svd(crv, 'econ');
            s = diag(s);
            r = my_chop2(s, tol*norm(s)/sqrt(d));
            v = conj(v(:,1:r));
            u = u(:,1:r)*diag(s(1:r));
            
            % Prepare enrichment, if needed
            if (kickrank>0)
                crv = u*v.';
                crv = reshape(crv, b, rv(i)*n(i)*rv(i+1));
                crv = crv(b,:).';
                % For updating z
                crz = bfun3(phizax{i}, A{i}, phizax{i+1}, crv);
                % For adding into solution
                crs = bfun3(phizax{i}, A{i}, phia{i+1}, crv);
                crs = reshape(crs, rz(i), n(i)*rv(i+1));
                v = [v,crs.'];
                [v,R]=qr(v, 0);
                u = [u, zeros(b*rv(i), rz(i))];
                u = u*R.';
                r = size(v, 2);
            end;
            cr2 = reshape(V{i-1}, rv(i-1)*n(i-1), rv(i));
            u = reshape(u, b, rv(i)*r);
            u = u.';
            u = reshape(u, rv(i), r*b);
            cr2 = cr2*u;
            V{i-1} = reshape(cr2, rv(i-1), n(i-1), r, b);
            V{i} = reshape(v.', r, n(i), rv(i+1));
           
            rv(i) = r;
            
            phia{i} = compute_next_Phi(phia{i+1}, V{i}, A{i}, V{i}, 'rl');
            phix{i} = compute_next_Phi(phix{i+1}, V{i}, [], x{i}, 'rl');
            
            if (kickrank>0)
                crz = reshape(crz, rz(i), n(i)*rz(i+1));
                [crz,R]=qr(crz.', 0);
                cr2 = reshape(z{i-1}, rz(i-1)*n(i-1), rz(i));
                cr2 = cr2*R.';
                rz(i) = size(crz, 2);
                z{i} = reshape(crz.', rz(i), n(i), rz(i+1));
                z{i-1} = reshape(cr2, rz(i-1), n(i-1), rz(i));
                
                phizax{i} = compute_next_Phi(phizax{i+1}, z{i}, A{i}, V{i}, 'rl');
            end;            
        end;
    end;
    
    if (verb>1)
        fprintf('amr-kryl: swp=[%d,%d], dx=%3.3e, r=%d, b=%d\n', swp, i, dx, r, b);
    end;
    
    % Stopping or reversing
    if ((dir>0)&&(i==d))||((dir<0)&&(i==1))
        if (verb>0)
            fprintf('amr-kryl: swp=%d{%d}, max_dx=%3.3e, max_r=%d, b=%d\n', swp, (1-dir)/2, max_dx, max(rv), b);
        end;
        if (max_dx<tol)&&(dir>0) % We must return V with b in V{d}
            break;
        else
            % We are at the terminal block
            % Here, we may allocate the space for a new krylov vector
            b = b_new;
            b_new = 0;
            crv = crv(:,1:b);
            if ((b+1)<=Mmax)
                crv = [crv, zeros(rv(i)*n(i)*rv(i+1), 1)];
                b = b+1;
            end;
            V{i} = reshape(crv, rv(i), n(i), rv(i+1), b);
            if (dir>0); swp = swp+1; end;
        end;
        dir = -dir;
        max_dx = 0;
    else
        i = i+dir;
    end;
end;
V{d} = reshape(crv, rv(d), n(d), b);
V = cell2core(tt_tensor, V);
if (abs(y(b))/abs(y(1))>tol)
    fprintf('!!! Krylov sequence has not converged (ymin/ymax=%3.1e).\n!!! Increase M or decrease |f(A)|\n', abs(y(b))/abs(y(1)));
end;
end


% new
function [Phi] = compute_next_Phi(Phi_prev, x, A, y, direction)
% Performs the recurrent Phi (or Psi) matrix computation
% Phi = Phi_prev * (x'Ay).
% If direction is 'lr', computes Psi
% if direction is 'rl', computes Phi
% A can be empty, then only x'y is computed.

% Phi1: rx1, ry1, ra1, or rx1, ry1
% Phi2: ry2, ra2, rx2, or ry2, rx2


rx1 = size(x,1); n = size(x,2); rx2 = size(x,3);
ry1 = size(y,1); m = size(y,2); ry2 = size(y,3);
if (~isempty(A))
    ra1 = size(A,1); ra2 = size(A,4);
else
    ra1 = 1; ra2 = 1;
end;

if (strcmp(direction, 'lr'))
    %lr: Phi1
    x = reshape(x, rx1, n*rx2);
    Phi = reshape(Phi_prev, rx1, ry1*ra1);
    Phi = x'*Phi;
    if (~isempty(A))
        Phi = reshape(Phi, n*rx2*ry1, ra1);
        Phi = Phi.';
        Phi = reshape(Phi, ra1*n, rx2*ry1);
        A = reshape(A, ra1*n, m*ra2);
        Phi = A.'*Phi;
        Phi = reshape(Phi, m, ra2*rx2*ry1);
    else
        Phi = reshape(Phi, n, rx2*ry1);
    end;
    Phi = Phi.';
    Phi = reshape(Phi, ra2*rx2, ry1*m);
    
    y = reshape(y, ry1*m, ry2);
    Phi = Phi*y;
    if (~isempty(A))
        Phi = reshape(Phi, ra2, rx2*ry2);
        Phi = Phi.';
    end;
    Phi = reshape(Phi, rx2, ry2, ra2);
else
    %rl: Phi2
    y = reshape(y, ry1*m, ry2);
    Phi = reshape(Phi_prev, ry2, ra2*rx2);
    Phi = y*Phi;
    if (~isempty(A))
        Phi = reshape(Phi, ry1, m*ra2*rx2);
        Phi = Phi.';
        Phi = reshape(Phi, m*ra2, rx2*ry1);
        A = reshape(A, ra1*n, m*ra2);
        Phi = A*Phi;
        Phi = reshape(Phi, ra1*n*rx2, ry1);
        Phi = Phi.';
    end;
    
    Phi = reshape(Phi, ry1*ra1, n*rx2);
    x = reshape(x, rx1, n*rx2);
    Phi = Phi*x';
    if (~isempty(A))
        Phi = reshape(Phi, ry1, ra1, rx1);
    else
        Phi = reshape(Phi, ry1, rx1);
    end;
end;


end


% new
function [y]=bfun3(Phi1, A, Phi2, x)
% Phi1: ry1, rx1, ra1
ry1 = size(Phi1,1);
rx1 = size(Phi1,2);
ra1 = size(Phi1,3);
% Phi2: rx2, ra2, ry2
ry2 = size(Phi2,3);
rx2 = size(Phi2,1);
ra2 = size(Phi2,2);
b = size(x,2);

n = size(A,2);
m = size(A,3);

y = reshape(x.', b*rx1*m, rx2);
Phi2 = reshape(Phi2, rx2, ra2*ry2);
y = y*Phi2;
y = reshape(y, b*rx1, m*ra2*ry2);
y = y.';
y = reshape(y, m*ra2, ry2*b*rx1);
A = reshape(A, ra1*n, m*ra2);
y = A*y;
y = reshape(y, ra1*n*ry2*b, rx1);
y = y.';
y = reshape(y, rx1*ra1, n*ry2*b);
Phi1 = reshape(Phi1, ry1, rx1*ra1);
y = Phi1*y;
y = reshape(y, ry1*n*ry2, b);
end