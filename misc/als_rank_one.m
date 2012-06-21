function [y] = als_rank_one(x, varargin)
% Computes rank-1 approximation to Ax-f via the ALS

d = x.d;
n = x.n;
A = [];
f = [];
tol = 1e-2;
maxit = 100;
verb = 0;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'a'
            A=varargin{i+1};
        case 'f'
            f=varargin{i+1};
        case 'tol'
            tol=varargin{i+1};
        case 'maxit'
            maxit=varargin{i+1};            
        case 'verb'
            verb=varargin{i+1};            
            
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if (isempty(A))
    A = tt_eye(n);
end;
if (isempty(f))
    f = tt_zeros(n);
end;

rf = f.r;

% y = tt_ones(n);
y = round(A*x-f, tol, 1);
% y = tt_rand(n, d, 1);
% y = y.*y;
% for i=1:d
%     y{i} = y{i} + 1;
% end;

cry = core2cell(y);
cra = core2cell(A);
crf = core2cell(f);
crx = core2cell(x);


phia = cell(d+1,1); phia{1}=1; phia{d+1}=1;
phif = cell(d+1,1); phif{1}=1; phif{d+1}=1;

for i=d:-1:2
    cry{i} = cry{i}./norm(cry{i}, 'fro');
    phia{i} = compute_next_Phi(phia{i+1}, cry{i}, cra{i}, crx{i}, 'rl');
    phif{i} = compute_next_Phi(phif{i+1}, cry{i}, [], crf{i}, 'rl');
end;

swp = 1;
i = 1;
dir = 1;
dx_max = 0;
while (swp<=maxit)
    Ax_loc = bfun3(phia{i}, cra{i}, phia{i+1}, crx{i});
    f_loc = phif{i}*reshape(crf{i}, rf(i), n(i)*rf(i+1));
    f_loc = reshape(f_loc, n(i), rf(i+1));
    f_loc = f_loc*(phif{i+1}.');
    
    newy = Ax_loc - f_loc;
    newy = reshape(newy, 1, n(i));
    dx = norm(newy-cry{i})/norm(newy);
    dx_max = max(dx, dx_max);
   
    nrm = norm(newy);
    newy = newy/norm(newy);
       
    if (dir==1)&&(i<d)
        cry{i+1} = cry{i+1}*nrm;
        phia{i+1} = compute_next_Phi(phia{i}, newy, cra{i}, crx{i}, 'lr');
        phif{i+1} = compute_next_Phi(phif{i}, newy, [], crf{i}, 'lr');
    elseif (dir==-1)&&(i>1)
        cry{i-1} = cry{i-1}*nrm;
        phia{i} = compute_next_Phi(phia{i+1}, newy, cra{i}, crx{i}, 'rl');
        phif{i} = compute_next_Phi(phif{i+1}, newy, [], crf{i}, 'rl');        
    else
        newy = newy*nrm;
    end;
    
    cry{i} = newy;
    
    i = i+dir;
    if ((i==d+1)&&(dir==1))||((i==0)&&(dir==-1))
        dir = -dir;
        
        if (verb>0)
            fprintf('--als_rank_one-- swp %d, dx_max: %3.3e\n', swp, dx_max);
        end;
        
        if (dx_max<tol)
            break;
        end;
        swp = swp+1;
        if (swp>maxit)&&(dx_max>tol)
            fprintf('--warn-- als_rank_one did not converge in %d iters\n', maxit);
        end;
        dx_max = 0;
        i = i+dir;
    end;
end;

y = cell2core(y, cry);
end


function [y]=bfun3(Phi1,B1,Phi2, x)
% Computes (Phi1 * B1 * Phi2)*x
% Phi1 is of sizes ry1, rB1, rx1
% B1 is of sizes rB1, k1, m1, rB2
% Phi2 is of sizes ry2, rB2, rx2
ry1 = size(Phi1,1); ry2 = size(Phi2,1);
rx1 = size(Phi1,3); rx2 = size(Phi2,3);
rb1=size(B1,1); rb2=size(B1,4); 
m1 = size(B1,3);
k1 = size(B1,2);

y = reshape(x, rx1, m1*rx2);
Phi1 = reshape(Phi1, ry1*rb1, rx1);
y = Phi1*y; % size ry1*rb1,m1*rx2 % cplx rb*rx^3*m^2
y = reshape(y, ry1, rb1*m1, rx2);
y = permute(y, [2, 1, 3]);
y = reshape(y, rb1*m1, ry1*rx2);
B1 = permute(B1, [2, 4, 1, 3]);
B1 = reshape(B1, k1*rb2, rb1*m1);
y = B1*y; % size k1*rb2, ry1*rx2 % cplx rb^2*rx^2*n^3
y = reshape(y, k1, rb2, ry1, rx2);
y = permute(y, [2, 4, 3, 1]);
y = reshape(y, rb2*rx2, ry1*k1);
Phi2 = reshape(Phi2, ry2, rb2*rx2);
y = Phi2*y; % size ry2, ry1*k1 % cplx rb*rx^3*n^2
y = y.';
y = reshape(y, ry1*k1*ry2, 1);
end


function [Phi] = compute_next_Phi(Phi_prev, x, A, y, direction)
% Performs the recurrent Phi (or Psi) matrix computation
% Phi = Phi_prev * (x'Ay).
% If direction is 'lr', computes Psi
% if direction is 'rl', computes Phi
% A can be empty, then only x'y is computed.

if (strcmp(direction, 'rl'))
  % Revert ranks to perform the right-to-left recursion
  x = permute(x, [3, 2, 1]);
  y = permute(y, [3, 2, 1]);
  if (~isempty(A))
    A = permute(A, [4, 2, 3, 1]);
  end
end

rx1 = size(x,1); n = size(x,2); rx2 = size(x,3);
ry1 = size(y,1); m = size(y,2); ry2 = size(y,3);
if (~isempty(A))
  ra1 = size(A,1); ra2 = size(A,4);
else
  ra1 = 1; ra2 = 1;
end

Phi = reshape(Phi_prev, [rx1*ra1, ry1]);
y = reshape(y, [ry1, m*ry2]);
Phi = Phi*y;	% complexity §\mcommentfont$\mathcal{O}(n  r_x r_A r_y^2)$§
Phi = reshape(Phi, [rx1, ra1, m, ry2]);
Phi = permute(Phi, [2, 3, 1, 4]);
if (~isempty(A))
  Phi = reshape(Phi, [ra1*m, rx1*ry2]);
  A = permute(A, [4, 2, 1, 3]);
  A = reshape(A, [ra2*n, ra1*m]);
  Phi = A*Phi;	% complexity §\mcommentfont$\mathcal{O}(n^2  r_x r_A^2 r_y)$§
  Phi = reshape(Phi, [ra2, n, rx1, ry2]);
end
Phi = permute(Phi, [3, 2, 1, 4]);
Phi = reshape(Phi, [rx1*n, ra2*ry2]);
x = reshape(x, [rx1*n, rx2]);
Phi = (x')*Phi;	% complexity §\mcommentfont$\mathcal{O}(n  r_x^2 r_A r_y)$§
if (~isempty(A))
  Phi = reshape(Phi, [rx2, ra2, ry2]);
end
end
