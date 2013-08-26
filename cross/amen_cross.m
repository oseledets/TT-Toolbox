function [y]=amen_cross(n, fun, tol, varargin)
% Block cross with error-based enrichment.
% function [y]=amen_cross(n, fun, tol, varargin)
%
% fun = @(i)fun(i)
% i comes as an array of sizes M x d,
% fun should return the block M x b (vectorized variant).
% M=1 if vec=false, hence the return may be either 1 x b or b x 1.

vars = varargin;

y = [];
nswp = 20;
kickrank = 4;
kickrank2 = 0;
tol_exit = tol;
verb = 1;
vec = false;

i = 1;
while (i<length(vars))
    switch lower(vars{i})
        case 'y0'
            y=vars{i+1};
        case 'nswp'
            nswp=vars{i+1};
        case 'kickrank'
            kickrank=vars{i+1};
        case 'kickrank2'
            kickrank2=vars{i+1};            
        case 'tol_exit'
            tol_exit=vars{i+1};  
        case 'verb'
            verb = vars{i+1};      
        case 'vec'
            vec = vars{i+1};               
        otherwise
            warning('Option %s was not recognized', vars{i});
    end;
    i=i+2;
end;

d = numel(n);
if (isempty(y))
    y = tt_rand(n, d, 2);
end;
ry = y.r;
y = core2cell(y);

if (kickrank>0)
    z = tt_rand(n,d,kickrank);
    rz = z.r;
    z = core2cell(z);
    % Interp matrices for z and y->z    
    phizz = cell(d+1,1);
    phizz{1} = 1;
    phizz{d+1} = 1;
    phizy = cell(d+1,1);
    phizy{1} = 1;
    phizy{d+1} = 1;
end;

% Interp matrices for y
phiyy = cell(d+1,1);
phiyy{1} = 1;
phiyy{d+1}=1;

Jy = cell(d+1,1);
Jz = cell(d+1,1);

% Initial orthog and indices: maxvol over y, z
for i=d:-1:2
    cry = reshape(y{i}, ry(i), n(i)*ry(i+1));
    [cry,rv]=qr(cry.', 0);
    cr2 = reshape(y{i-1}, ry(i-1)*n(i-1), ry(i));
    cr2 = cr2*rv.';
    ry(i) = size(cry,2);
    y{i} = reshape(cry.', ry(i), n(i), ry(i+1));
    y{i-1} = reshape(cr2, ry(i-1), n(i-1), ry(i));
    
    cry = reshape(y{i}, ry(i)*n(i), ry(i+1));
    cry = cry*phiyy{i+1};
    cry = reshape(cry, ry(i), n(i)*ry(i+1));
    ind = maxvol2(cry.');
    phiyy{i} = cry(:,ind);
    Jy{i} = [kron(ones(ry(i+1),1), (1:n(i))'), kron(Jy{i+1}, ones(n(i),1))]; % n*r2, d+1-i
    Jy{i} = Jy{i}(ind,:);
    
    if (kickrank>0)
        % The same for residual
        crz = reshape(z{i}, rz(i), n(i)*rz(i+1));
        [crz,rv]=qr(crz.', 0);
        cr2 = reshape(z{i-1}, rz(i-1)*n(i-1), rz(i));
        cr2 = cr2*rv.';
        rz(i) = size(crz,2);
        z{i} = reshape(crz.', rz(i), n(i), rz(i+1));
        z{i-1} = reshape(cr2, rz(i-1), n(i-1), rz(i));
        
        crz = reshape(z{i}, rz(i)*n(i), rz(i+1));
        crz = crz*phizz{i+1};
        crz = reshape(crz, rz(i), n(i)*rz(i+1));
        ind = maxvol2(crz.');
        phizz{i} = crz(:,ind);
        Jz{i} = [kron(ones(rz(i+1),1), (1:n(i))'), kron(Jz{i+1}, ones(n(i),1))]; % n*r2, d+1-i
        Jz{i} = Jz{i}(ind,:);
        % Extract z indices from y.
        cry = reshape(y{i}, ry(i)*n(i), ry(i+1));
        cry = cry*phizy{i+1};
        cry = reshape(cry, ry(i), n(i)*rz(i+1));
        phizy{i} = cry(:, ind);
    end;
end;

evalcnt = 0;
b = []; % A block size will be stored here
max_dx = 0;
swp = 1;
dir = 1;
i = 1;
while (swp<=nswp)
    % Generate a full index = [Jl, (1:n)', Jr]
    % Mind big-endian: Jyl to lowest index           n(i) to the middle one,                      Jyr to the senior index
    J = [kron(ones(n(i)*ry(i+1),1), Jy{i}), kron(kron(ones(ry(i+1),1), (1:n(i))'), ones(ry(i),1)), kron(Jy{i+1}, ones(ry(i)*n(i),1))];
    evalcnt = evalcnt + size(J,1);
    if (isempty(b))
        % Take one sample to measure its size =)
        cry1 = fun(J(1,:));
        b = numel(cry1);
        % Init y{1} by zeros
        y{1} = zeros(ry(1), n(1), ry(2), b);
        evalcnt = evalcnt+1;
    end;
    % Actually compute the new core
    if (vec)
        cry = fun(J);
    else
        % User function is not vectorized - run a loop
        cry = zeros(ry(i)*n(i)*ry(i+1), b);
        for j=1:ry(i)*n(i)*ry(i+1)
            cry(j,:) = fun(J(j,:));
        end;
    end;
    
    % Apply interpolation matrices
    cry = reshape(cry, ry(i), n(i)*ry(i+1)*b);
    cry = phiyy{i} \ cry;
    cry = reshape(cry, ry(i)*n(i)*ry(i+1), b);
    cry = cry.';
    cry = reshape(cry, b*ry(i)*n(i), ry(i+1));
    cry = cry/phiyy{i+1};
    cry = reshape(cry, b, ry(i)*n(i)*ry(i+1));
    cry = cry.';
    
    y_prev = reshape(y{i}, ry(i)*n(i)*ry(i+1), b);
    dx = norm(cry-y_prev, 'fro')/norm(cry, 'fro');
    max_dx = max(max_dx,dx);
    
    % Truncation, etc
    if (dir>0)&&(i<d)
        cry = reshape(cry, ry(i)*n(i), ry(i+1)*b);
        [u,s,v]=svd(cry, 'econ');
        s = diag(s);
        r = my_chop2(s, norm(s)*tol/sqrt(d));
        u = u(:,1:r);
        v = diag(s(1:r))*v(:,1:r)';
        
        if (kickrank>0)
            % Project onto Ez
            crys = u*v;
            crys = reshape(crys, ry(i)*n(i)*ry(i+1), b);
            crys = crys.';
            crys = reshape(crys, b*ry(i)*n(i), ry(i+1));
            crys = crys*phizy{i+1};
            crys = reshape(crys, b, ry(i)*n(i)*rz(i+1));
            crys = crys.';
            cryz = reshape(crys, ry(i), n(i)*rz(i+1)*b);
            cryz = phizy{i}*cryz;
            cryz = reshape(cryz, rz(i)*n(i)*rz(i+1), b);
            
            % Compute Z
            J = [kron(ones(n(i)*rz(i+1),1), Jz{i}), kron(kron(ones(rz(i+1),1), (1:n(i))'), ones(rz(i),1)), kron(Jz{i+1}, ones(rz(i)*n(i),1))];
            evalcnt = evalcnt + size(J,1);
            if (vec)
                crz = fun(J);
            else
                % User function is not vectorized - run a loop
                crz = zeros(rz(i)*n(i)*rz(i+1), b);
                for j=1:rz(i)*n(i)*rz(i+1)
                    crz(j,:) = fun(J(j,:));
                end;
            end;    
            crz = crz - cryz;            
            % Apply interpolation matrices
            crz = reshape(crz, rz(i), n(i)*rz(i+1)*b);
            crz = phizz{i} \ crz;
            crz = reshape(crz, rz(i)*n(i)*rz(i+1), b);
            crz = crz.';
            crz = reshape(crz, b*rz(i)*n(i), rz(i+1));
            crz = crz/phizz{i+1};
            crz = reshape(crz, b, rz(i)*n(i)*rz(i+1));
            crz = crz.';
            crz = reshape(crz, rz(i)*n(i), rz(i+1)*b);
            [crz,sz,vz]=svd(crz, 'econ');
            crz = crz(:,1:min(size(crz,2),kickrank));
            if (kickrank2>0)
                crz = [crz, randn(rz(i)*n(i), kickrank2)];
                [crz,sz]=qr(crz,0);
            end;
            
            % Compute S
            J = [kron(ones(n(i)*rz(i+1),1), Jy{i}), kron(kron(ones(rz(i+1),1), (1:n(i))'), ones(ry(i),1)), kron(Jz{i+1}, ones(ry(i)*n(i),1))];
            evalcnt = evalcnt + size(J,1);
            if (vec)
                crs = fun(J);
            else
                % User function is not vectorized - run a loop
                crs = zeros(ry(i)*n(i)*rz(i+1), b);
                for j=1:ry(i)*n(i)*rz(i+1)
                    crs(j,:) = fun(J(j,:));
                end;
            end;    
            % Apply interpolation matrices - here crys is already in the
            % core variables from the left
            crs = reshape(crs, ry(i), n(i)*rz(i+1)*b);
            crs = phiyy{i} \ crs;
            crs = reshape(crs, ry(i)*n(i)*rz(i+1), b);
            
            crs = crs - crys;
            
            crs = crs.';
            crs = reshape(crs, b*ry(i)*n(i), rz(i+1));
            crs = crs/phizz{i+1};
            crs = reshape(crs, b, ry(i)*n(i)*rz(i+1));
            crs = crs.';
            crs = reshape(crs, ry(i)*n(i), rz(i+1)*b);
            [crs,sz,vz]=svd(crs, 'econ');
            crs = crs(:,1:min(size(crs,2),kickrank));
            
            % Enrichment itself
            [u,rv]=qr([u,crs], 0);
            v = [v; zeros(size(crs,2), ry(i+1)*b)];
            v = rv*v;
        end;
        
        r = size(u,2);
        v = reshape(v, r*ry(i+1), b);
        v = v.';
        v = reshape(v, b*r, ry(i+1));
        cr2 = reshape(y{i+1}, ry(i+1), n(i+1)*ry(i+2));
        cr2 = v*cr2;
        cr2 = reshape(cr2, b, r*n(i+1)*ry(i+2));
        cr2 = cr2.';
        ry(i+1) = r;        
        y{i} = reshape(u, ry(i), n(i), ry(i+1));
        y{i+1} = reshape(cr2, ry(i+1), n(i+1), ry(i+2), b);
        
        % Maxvols and phis
        u = reshape(u, ry(i), n(i)*ry(i+1));
        u = phiyy{i}*u;
        u = reshape(u, ry(i)*n(i), ry(i+1));
        ind = maxvol2(u);
        phiyy{i+1} = u(ind,:);
        Jy{i+1} = [kron(ones(n(i),1), Jy{i}), kron((1:n(i))', ones(ry(i),1))];
        Jy{i+1} = Jy{i+1}(ind,:);
        if (kickrank>0)
            rz(i+1) = size(crz,2);
            z{i} = crz;
            crz = reshape(crz, rz(i), n(i)*rz(i+1));
            crz = phizz{i}*crz;
            crz = reshape(crz, rz(i)*n(i), rz(i+1));
            ind = maxvol2(crz);
            phizz{i+1} = crz(ind,:);
            Jz{i+1} = [kron(ones(n(i),1), Jz{i}), kron((1:n(i))', ones(rz(i),1))];
            Jz{i+1} = Jz{i+1}(ind,:);
            phizy{i+1} = reshape(y{i}, ry(i), n(i)*ry(i+1));
            phizy{i+1} = phizy{i}*phizy{i+1};
            phizy{i+1} = reshape(phizy{i+1}, rz(i)*n(i), ry(i+1));
            phizy{i+1} = phizy{i+1}(ind,:);
        end;
    elseif (dir<0)&&(i>1)
        cry = cry.';
        cry = reshape(cry, b*ry(i), n(i)*ry(i+1));
        [u,s,v]=svd(cry, 'econ');
        s = diag(s);
        r = my_chop2(s, norm(s)*tol/sqrt(d));
        u = u(:,1:r)*diag(s(1:r));
        v = v(:,1:r)';
        
        if (kickrank>0)
            % Project onto Ez
            crys = u*v;
            crys = reshape(crys, b, ry(i)*n(i)*ry(i+1));
            crys = crys.';
            crys = reshape(crys, ry(i), n(i)*ry(i+1)*b);
            crys = phizy{i}*crys;
            crys = reshape(crys, rz(i)*n(i)*ry(i+1), b);
            crys = crys.';
            cryz = reshape(crys, b*rz(i)*n(i), ry(i+1));
            cryz = cryz*phizy{i+1};
            cryz = reshape(cryz, b, rz(i)*n(i)*rz(i+1));
            
            % Compute Z
            J = [kron(ones(n(i)*rz(i+1),1), Jz{i}), kron(kron(ones(rz(i+1),1), (1:n(i))'), ones(rz(i),1)), kron(Jz{i+1}, ones(rz(i)*n(i),1))];
            evalcnt = evalcnt + size(J,1);
            if (vec)
                crz = fun(J);
                crz = crz.';
            else
                % User function is not vectorized - run a loop
                crz = zeros(b, rz(i)*n(i)*rz(i+1));
                for j=1:rz(i)*n(i)*rz(i+1)
                    crz(:,j) = fun(J(j,:));
                end;
            end;    
            crz = crz - cryz;            
            % Apply interpolation matrices
            crz = reshape(crz, b*rz(i)*n(i), rz(i+1));
            crz = crz/phizz{i+1};            
            crz = reshape(crz, b, rz(i)*n(i)*rz(i+1));
            crz = crz.';
            crz = reshape(crz, rz(i), n(i)*rz(i+1)*b);
            crz = phizz{i} \ crz;            
            crz = reshape(crz, rz(i)*n(i)*rz(i+1), b);
            crz = crz.';
            crz = reshape(crz, b*rz(i), n(i)*rz(i+1));
            [vz,sz,crz]=svd(crz, 'econ');
            crz = conj(crz(:,1:min(size(crz,2),kickrank)));
            if (kickrank2>0)
                crz = [crz, randn(n(i)*rz(i+1), kickrank2)];
                [crz,sz]=qr(crz,0);
            end;            
            crz = crz.';
            
            % Compute S
            J = [kron(ones(n(i)*ry(i+1),1), Jz{i}), kron(kron(ones(ry(i+1),1), (1:n(i))'), ones(rz(i),1)), kron(Jy{i+1}, ones(rz(i)*n(i),1))];
            evalcnt = evalcnt + size(J,1);
            if (vec)
                crs = fun(J);
                crs = crs.';
            else
                % User function is not vectorized - run a loop
                crs = zeros(b, rz(i)*n(i)*ry(i+1));
                for j=1:rz(i)*n(i)*ry(i+1)
                    crs(:,j) = fun(J(j,:));
                end;
            end;    
            % Apply interpolation matrices
            crs = reshape(crs, b*rz(i)*n(i), ry(i+1));
            crs = crs/phiyy{i+1};            
            crs = reshape(crs, b, rz(i)*n(i)*ry(i+1));
            
            crs = crs - crys;

            crs = crs.';
            crs = reshape(crs, rz(i), n(i)*ry(i+1)*b);
            crs = phizz{i} \ crs;            
            crs = reshape(crs, rz(i)*n(i)*ry(i+1), b);
            crs = crs.';
            crs = reshape(crs, b*rz(i), n(i)*ry(i+1));
            [vz,sz,crs]=svd(crs, 'econ');
            crs = crs(:,1:min(size(crs,2),kickrank))';
            
            % Enrichment itself
            [v,rv]=qr([v;crs].', 0);
            u = [u, zeros(b*ry(i), size(crs,1))];
            u = u*rv.';
            v = v.';
        end;
        
        r = size(v,1);
        u = reshape(u, b, ry(i)*r);
        u = u.';
        u = reshape(u, ry(i), r*b);
        cr2 = reshape(y{i-1}, ry(i-1)*n(i-1), ry(i));
        cr2 = cr2*u;
        ry(i) = r;        
        y{i} = reshape(v, ry(i), n(i), ry(i+1));
        y{i-1} = reshape(cr2, ry(i-1), n(i-1), ry(i), b);
        
        % Maxvols and phis
        v = reshape(v, ry(i)*n(i), ry(i+1));
        v = v*phiyy{i+1};
        v = reshape(v, ry(i), n(i)*ry(i+1));
        ind = maxvol2(v.');
        phiyy{i} = v(:,ind);
        Jy{i} = [kron(ones(ry(i+1),1), (1:n(i))'), kron(Jy{i+1}, ones(n(i),1))]; % n*r2, d+1-i
        Jy{i} = Jy{i}(ind,:);        
        if (kickrank>0)
            rz(i) = size(crz,1);
            z{i} = crz;
            crz = reshape(crz, rz(i)*n(i), rz(i+1));
            crz = crz*phizz{i+1};
            crz = reshape(crz, rz(i), n(i)*rz(i+1));
            ind = maxvol2(crz.');
            phizz{i} = crz(:,ind);
            Jz{i} = [kron(ones(rz(i+1),1), (1:n(i))'), kron(Jz{i+1}, ones(n(i),1))]; % n*r2, d+1-i
            Jz{i} = Jz{i}(ind,:);
            % Extract indices from y.
            cry = reshape(y{i}, ry(i)*n(i), ry(i+1));
            cry = cry*phizy{i+1};
            cry = reshape(cry, ry(i), n(i)*rz(i+1));
            phizy{i} = cry(:, ind);
        end;
    else
        y{i} = reshape(cry, ry(i), n(i), ry(i+1), b);
    end;
    
    i = i+dir;
    
    % Check the convergence, restart
    if ((i==d+1)||(i==0))
        if (verb>0)
            fprintf('=amen_cross= swp=%d, max_dx=%3.3e, max_rank=%d, cum#evals=%d\n', swp, max_dx, max(ry), evalcnt);
        end;
        
        if (dir>0)
            if (max_dx<tol_exit)
                break;
            end;
        end;
        
        max_dx = 0;
        dir = -dir;
        i = i+dir;
        swp = swp+1;
    end;
end;

y{d} = reshape(y{d}, ry(d), n(d), b);
y = cell2core(tt_tensor,y);
end
