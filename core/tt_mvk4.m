function [y]=tt_mvk4(A, x, tol, varargin)

% step_dpow = 0.1; % stepsize for d-power in truncations
% min_dpow = 1; % Minimal d-power for truncation
%
% bot_conv = 0.1; % bottom convergence factor - if better, we can decrease dpow and drank
% top_conv = 0.99; % top convergence factor - if worse, we have to increase dpow and drank


nswp=10;

als_tol_low = 5;
als_iters = 4;

rmax=1000;

verb=1;
kickrank = 10;
 kicktype = 'mr';
% kicktype = 'rand';

% pcatype = 'svd';
pcatype = 'uchol';

y=[];
block_order = [];

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
            rmax=varargin{i+1};
        case 'kicktype'
            kicktype=lower(varargin{i+1});
        case 'pcatype'
            pcatype=lower(varargin{i+1});            
        case 'y0'
            y=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case 'step_dpow'
            step_dpow=varargin{i+1};
        case 'min_dpow'
            min_dpow=varargin{i+1};
        case 'bot_conv'
            bot_conv=varargin{i+1};
        case 'top_conv'
            top_conv=varargin{i+1};
        case 'block_order'
            block_order=varargin{i+1};
        case 'als_tol_low'
            als_tol_low=varargin{i+1};
        case 'als_iters'
            als_iters=varargin{i+1};

        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

rx = x.r;
d = x.d;
tailranks = [false,false];
if (rx(1)~=1)
    e1 = tt_tensor;
    e1.d=1;
    e1.n = [rx(1)];
    e1.r = [1; rx(1)];
    e1.ps = [1; rx(1)^2+1];
    e1.core = reshape(eye(rx(1)), rx(1)^2, 1);
    x = kron(e1, x);
    
    A = kron(tt_matrix(eye(rx(1))), A);
    tailranks(1)=true;
end;
if (rx(d+1)~=1)
    e1 = tt_tensor;
    e1.d=1;
    e1.n = [rx(d+1)];
    e1.r = [rx(d+1); 1];
    e1.ps = [1; rx(d+1)^2+1];
    e1.core = reshape(eye(rx(d+1)), rx(d+1)^2, 1);
    x = kron(x, e1);
    
    A = kron(A, tt_matrix(eye(rx(d+1))));
    tailranks(2)=true;
end;

d = x.d;
rx = x.r;
ra = A.r;
n = A.n;
m = A.m;
if (isempty(y))
%     y = tt_rand(n, A.d, kickrank);
    y = tt_ones(n, A.d);
end;

if (isempty(block_order))
    block_order = [+(d), -(d)];
end;

ry = y.r;

cry = core2cell(y);
crA = core2cell(A);
crx = core2cell(x);

phia = cell(d+1,1); phia{1}=1; phia{d+1}=1;


% Orthogonalization
for i=d:-1:2
    cr = cry{i};
    cr = reshape(cr, ry(i), n(i)*ry(i+1));
    [cr, rv]=qr(cr.', 0);
    cr2 = cry{i-1};
    cr2 = reshape(cr2, ry(i-1)*n(i-1), ry(i));
    cr2 = cr2*(rv.');
    ry(i) = size(cr, 2);
    cr = reshape(cr.', ry(i), n(i), ry(i+1));
%     psy(i) = psy(i+1)-ry(i)*n(i)*ry(i+1);
%     cry(psy(i):(psy(i+1)-1))=cr(:);
%     psyo(i) = psyo(i-1)+ry(i-1)*n(i-1)*ry(i);
%     cry(psyo(i-1):(psyo(i)-1)) = cr2(:);
    cry{i-1} = reshape(cr2, ry(i-1), n(i-1), ry(i));
    cry{i} = cr;

    phia{i} = compute_next_Phi(phia{i+1}, cr, crA{i}, crx{i}, 'rl');
%     phia{i} = compute_next_Phi(phia{i+1}, cr, reshape(cra(psa(i):(psa(i+1)-1)), ra(i), n(i), n(i), ra(i+1)), reshape(crx(psx(i):(psx(i+1)-1)), rx(i), n(i), rx(i+1)), 'rl');

end;
% cry = cry(psy(1):(psy(d+1)-1));

last_sweep = false;
swp = 1;
i = 1;

% dy_old = ones(d,1);
dy = zeros(d,1);
max_dy = 0;
% For extra-rank addition
% dpows = ones(d,1)*min_dpow;
% dranks = zeros(d,1);

cur_order = block_order;
order_index = 1;
dir = sign(cur_order(order_index));

% DMRG sweeps
while (swp<=nswp)

    % Extract elements - matrix
    Phi1 = phia{i}; Phi2 = phia{i+1};
%     A1 = reshape(cra(psa(i):(psa(i+1)-1)), ra(i), n(i), n(i), ra(i+1));
%     x1 = reshape(crx(psx(i):(psx(i+1)-1)), rx(i), n(i), rx(i+1));
    x1 = crx{i};
    A1 = crA{i};

    newy = reshape(Phi1, ry(i)*ra(i), rx(i));
    newy = newy*reshape(x1, rx(i), m(i)*rx(i+1));
    newy = reshape(newy, ry(i), ra(i)*m(i), rx(i+1));
    newy = permute(newy, [2, 1, 3]);
    newy = reshape(newy, ra(i)*m(i), ry(i)*rx(i+1));
    A1 = permute(A1, [2, 4, 1, 3]);
    A1 = reshape(A1, n(i)*ra(i+1), ra(i)*m(i));
    newy = A1*newy;
    newy = reshape(newy, n(i), ra(i+1), ry(i), rx(i+1));
    newy = permute(newy, [3, 1, 2, 4]);
    newy = reshape(newy, ry(i)*n(i), ra(i+1)*rx(i+1));
    if (dir>0)&&(strcmp(kicktype, 'mr'))
        newy_save = newy;
    end;
    newy = newy*(reshape(Phi2, ry(i+1), ra(i+1)*rx(i+1)).');
    newy = reshape(newy, ry(i)*n(i)*ry(i+1), 1);

    oldy = reshape(cry{i}, ry(i)*n(i)*ry(i+1), 1);

    dy(i) = norm(newy-oldy)/norm(newy);
    max_dy = max(max_dy, dy(i));
%
%     % The new core does not converge - increase rank
%     if (dy(i)/dy_old(i)>top_conv)&&(dy(i)>tol)
%         dranks(i)=dranks(i)+1;
%         dpows(i)=dpows(i)+step_dpow;
%     end;
%     % The new core converges well - try to decrease rank
%     if (dy(i)/dy_old(i)<bot_conv)||(dy(i)<tol)
%         dranks(i)=max(dranks(i)-1, 0);
%         dpows(i)=max(dpows(i)-step_dpow, min_dpow);
%     end;
%
%     if (last_sweep)
%         dpows(i)=0.5;
%         dranks(i)=0;
%     end;


    % Truncation
    if (dir>0) % left-to-right
        newy = reshape(newy, ry(i)*n(i), ry(i+1));
    else
        newy = reshape(newy, ry(i), n(i)*ry(i+1));
    end;

    if (kickrank>=0)
    [u,s,v]=svd(newy, 'econ');
    s = diag(s);
    r = my_chop2(s, tol/sqrt(d)*norm(s));
%     % Artificial rank increasing
%     r = r+dranks(i);
    r = min(r, numel(s));
    r = min(r, rmax);
    else
        if (dir>0)
            [u,v]=qr(newy, 0);
            v=v';
            r = size(u,2);
            s = ones(r,1);
        else
            [v,u]=qr(newy.', 0);
            v=conj(v);
            u=u.';
            r = size(u,2);
            s = ones(r,1);
        end;
    end;

    if (verb>1)
	fprintf('=mvk4=   block %d{%d}, dy: %3.3e, r: %d\n', i, dir, dy(i), r);
%     elseif (verb>0)
%         fprintf('                                                                          \r')
%         fprintf('=mvk4=   block %d{%d}, dy: %3.3e, r: %d\r', i, dir, dy(i), r);
    end;

    if (dir>0)&&(i<d) % left-to-right, kickrank, etc
        u = u(:,1:r);
        v = v(:,1:r)*diag(s(1:r));

        % kick
        radd = 0; rv = 1;
        if (~last_sweep)&&(kickrank>0)
%             newy = reshape(Phi1, ry(i)*ra(i), rx(i));
%             newy = newy*reshape(crx{i}, rx(i), n(i)*rx(i+1));
%             newy = reshape(newy, ry(i), ra(i)*n(i), rx(i+1));
%             newy = permute(newy, [2, 1, 3]);
%             newy = reshape(newy, ra(i)*n(i), ry(i)*rx(i+1));
%             A1 = permute(crA{i}, [2, 4, 1, 3]);
%             A1 = reshape(A1, n(i)*ra(i+1), ra(i)*n(i));
%             newy = A1*newy;
%             newy = reshape(newy, n(i), ra(i+1), ry(i), rx(i+1));
%             newy = permute(newy, [3, 1, 2, 4]);
%             newy = reshape(newy, ry(i)*n(i), ra(i+1)*rx(i+1));
            if (strcmp(kicktype, 'mr'))
                leftresid = [reshape(u*v', ry(i)*n(i), ry(i+1)), -newy_save];
                if (strcmp(pcatype, 'svd'))
                    [uk,sk,vk]=svd(leftresid, 'econ');
                    uk = uk(:,1:min(kickrank, size(uk,2)));
                else
                    uk = uchol(leftresid.', kickrank*2);
                    uk = uk(:,size(uk,2):-1:max(size(uk,2)-kickrank+1,1));
                end;

%                 if (i+1<d)
%                     newy = reshape(phia{i+3}, ry(i+3)*ra(i+3), rx(i+3));
%                     newy = newy*reshape(crx{i+2}, rx(i+2)*n(i+2), rx(i+3)).';
%                     newy = reshape(newy, ry(i+3), ra(i+3), rx(i+2), n(i+2));
%                     newy = permute(newy, [2, 4, 1, 3]);
%                     newy = reshape(newy, ra(i+3)*n(i+2), ry(i+3)*rx(i+2));
%                     A1 = permute(crA{i+2}, [2, 1, 4, 3]);
%                     A1 = reshape(A1, n(i+2)*ra(i+2), ra(i+3)*n(i+2));
%                     newy = A1*newy;
%                     newy = reshape(newy, n(i+2), ra(i+2), ry(i+3), rx(i+2));
%                     newy = permute(newy, [1, 3, 2, 4]);
%                     newy = reshape(newy, n(i+2)*ry(i+3), ra(i+2)*rx(i+2));
%
%                     oldy = reshape(cry{i+2}, ry(i+2), n(i+2)*ry(i+3)).';
%                     rightresid = [oldy, -newy];
%                     vk = uchol(rightresid.', kickrank*2);
%                     vk = vk(:,size(vk,2):-1:max(size(vk,2)-kickrank+1,1));
%
%                     [oldy,rv]=qr([oldy,vk], 0);
%                     newy = reshape(cry{i+1}, ry(i+1)*n(i+1), ry(i+2));
%                     newy = [newy, zeros(ry(i+1)*n(i+1), size(vk,2))];
%                     newy = newy*(rv.');
%                     ry(i+2) = size(oldy,2);
%                     cry{i+1} = reshape(newy, ry(i+1), n(i+1), ry(i+2));
%                     cry{i+2} = permute(reshape(oldy, n(i+2), ry(i+3), ry(i+2)), [3, 1, 2]);
%                     phia{i+2} = compute_next_Phi(phia{i+3}, cry{i+2}, crA{i+2}, crx{i+2}, 'rl');
%                 end;
            else
                uk = rand(ry(i)*n(i), kickrank);
            end;
            [u,rv]=qr([u,uk], 0);
            radd = size(uk,2);
        end;
%         radd = size(u, 2)+size(uk,2)-r;
        v = [v, zeros(ry(i+1), radd)];
        v = rv*(v');
        r = size(u,2);
%         r = r+radd;
        cr2 = cry{i+1};
        cr2 = reshape(cr2, ry(i+1), n(i+1)*ry(i+2));
        v = v*cr2; % size r+radd, n2, r3

%         y{i} = reshape(u, ry(i), n(i), r);
%         y{i+1} = reshape(v, r, n(i+1), ry(i+2));
        ry(i+1) = r;

%         % Add residual component
%         rescomp = als_rank_one(x, 'A', A, 'f', y);
%         % We have to update ALL phis
%         Rl = [1,1];
%         for j=1:i
%             addbas = rescomp{i};
%             addbas = addbas/norm(addbas);
%             newbas = zeros(ry(i)+1, n(i), ry(i+1)+1);
%             newbas(1:ry(i), :, 1:ry(i+1))=y{i};
%             newbas(ry(i)+1, :, ry(i+1)+1)=addbas;
%             newbas = reshape(newbas, ry(i)+1, n(i)*(ry(i+1)+1));
%             newbas = Rl*newbas;
%             ry(i) = size(Rl, 1);
%             newbas = reshape(newbas, ry(i)*n(i), ry(i+1)+1);
%             % We have to orthogonalize it. approximately, if j<i
%             proj = newbas(:,ry(i+1)+1) - newbas(:,1:ry(i+1))*(newbas(:,1:ry(i+1))'*newbas(:,ry(i+1)+1));
%             if (i<j)&&(norm(proj)<tol/(d^dpows(i)))
%                 % Remove the new vector
%                 newbas2 = newbas(:, 1:ry(i+1));
%                 Rl = [eye(ry(i+1)); zeros(1,ry(i+1))];
%             else
%                 % Reort the new vector to the old one
%                 newbas2 = reort(newbas(:,1:ry(i+1)), newbas(:, ry(i+1)+1));
%                 Rl = newbas2'*newbas;
%             end;
%             radd = size(newbas2, 2)-ry(i+1);
%             y{i} = reshape(newbas2, ry(i), n(i), ry(i+1)+radd);
%             if (radd>0) % update phi
%                 addbas = reshape(newbas2(:, ry(i+1)+1), ry(i), n(i), 1);
%                 newphia = zeros(ry(i+1)+1, ra(i+1), rx(i+1));
%                 newphia(1:ry(i+1),:,:)=phia{i+1};
%                 newphia(ry(i+1)+1,:,:) = compute_next_Phi(phia{i}, addbas, A{i}, x{i}, 'lr');
%                 phia{i+1} = newphia;
%             end;
%         end;
%
%         Rr = [1;1];
%         for j=d:-1:i+2
%             addbas = rescomp{i};
%             addbas = addbas/norm(addbas);
%             newbas = zeros(ry(i)+1, n(i), ry(i+1)+1);
%             newbas(1:ry(i), :, 1:ry(i+1))=y{i};
%             newbas(ry(i)+1, :, ry(i+1)+1)=addbas;
%             newbas = reshape(newbas, (ry(i)+1)*n(i), (ry(i+1)+1));
%             newbas = newbas*Rr;
%             ry(i+1) = size(Rr, 2);
%             newbas = reshape(newbas, ry(i)+1, n(i)*ry(i+1));
%             newbas = newbas.';
%             % We have to orthogonalize it. approximately, if j<i
%             proj = newbas(:, ry(i)+1) - newbas(:,1:ry(i))*(newbas(:,1:ry(i))'*newbas(:,ry(i)+1));
%             if (i<j)&&(norm(proj)<tol/(d^dpows(i)))
%                 % Remove the new vector
%                 newbas2 = newbas(:, 1:ry(i));
%                 Rr = [eye(ry(i)), zeros(ry(i), 1)];
%             else
%                 % Reort the new vector to the old one
%                 newbas2 = reort(newbas(:,1:ry(i)), newbas(:, ry(i)+1));
%                 Rr = newbas2'*newbas;
%                 Rr = Rr.';
%             end;
%             radd = size(newbas2, 2)-ry(i);
%             y{i} = reshape(newbas2.', ry(i)+radd, n(i), ry(i+1));
%             if (radd>0) % update phi
%                 addbas = reshape(newbas2(:, ry(i)+1).', 1, n(i), ry(i+1));
%                 newphia = zeros(ry(i)+1, ra(i), rx(i));
%                 newphia(1:ry(i),:,:)=phia{i};
%                 newphia(ry(i)+1,:,:) = compute_next_Phi(phia{i+1}, addbas, A{i}, x{i}, 'lr');
%                 phia{i} = newphia;
%             end;
%         end;

        u = reshape(u, ry(i), n(i), r);
        v = reshape(v, r, n(i+1), ry(i+2));

        % Recompute phi. Left ones, so permute them appropriately
        phia{i+1} = compute_next_Phi(phia{i}, u, crA{i}, crx{i}, 'lr');

        % Stuff back

        cry{i} = u;
        cry{i+1} = v;
    elseif (dir<0)&&(i>1) % right-to-left
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        % kick
        radd = 0; rv = 1;
%          if (~last_sweep)&&(strcmp(kicktype, 'rand'))&&(kickrank>0)
%  %             newy = reshape(Phi2, ry(i+1)*ra(i+1), rx(i+1));
%  %             newy = newy*reshape(crx{i}, rx(i)*n(i), rx(i+1)).';
%  %             newy = reshape(newy, ry(i+1), ra(i+1), rx(i), n(i));
%  %             newy = permute(newy, [2, 4, 1, 3]);
%  %             newy = reshape(newy, ra(i+1)*n(i), ry(i+1)*rx(i));
%  %             A1 = permute(crA{i}, [1, 2, 4, 3]);
%  %             A1 = reshape(A1, ra(i)*n(i), ra(i+1)*n(i));
%  %             newy = A1*newy;
%  %             newy = reshape(newy, ra(i), n(i), ry(i+1), rx(i));
%  %             newy = permute(newy, [1, 4, 2, 3]);
%  %             newy = reshape(newy, ra(i)*rx(i), n(i)*ry(i+1));
%  %             resid = [reshape(u*v', ry(i), n(i)*ry(i+1)); -newy];
%  %             uk = uchol(resid, kickrank+1);
%  %             uk = uk(:,end:-1:max(end-kickrank+1,1));
%              uk = randn(n(i)*ry(i+1), kickrank);
%              [v,rv]=qr([v,uk], 0);
%              radd = size(uk,2);
%  %             v = reort(v, uk);
%          end;

%  %         radd = size(v, 2)+size(uk,2)-r;
        u = [u, zeros(ry(i), radd)];
        u = u*(rv.');
        cr2 = cry{i-1};
        cr2 = reshape(cr2, ry(i-1)*n(i-1), ry(i));
        u = cr2*u;

%         r = r+radd;
        r = size(v,2);

        u = reshape(u, ry(i-1), n(i-1), r);
        v = reshape(v.', r, n(i), ry(i+1));

        % Recompute phi. Here are right phis
        phia{i} = compute_next_Phi(phia{i+1}, v, crA{i}, crx{i}, 'rl');

        % Stuff back
        ry(i) = r;
        cry{i-1} = u;
        cry{i} = v;
    elseif ((dir>0)&&(i==d))||((dir<0)&&(i==1))
        % Just stuff back the last core
        newy = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
        newy = reshape(newy, ry(i), n(i), ry(i+1));
        cry{i} = newy;
    end;


    i = i+dir;

    % Reversing, residue check, etc
    cur_order(order_index) = cur_order(order_index) - dir;
    % New direction
    if (cur_order(order_index)==0)
        order_index = order_index+1;

        if (verb>0)
%  	    fprintf('                                                   \r')
            fprintf('=mvk4= sweep %d{%d}, max_dy: %3.3e, erank: %g\n', swp, order_index-1, max_dy, sqrt(ry(1:d)'*(n.*ry(2:d+1))/sum(n)));
        end;

        if (last_sweep)
            break;
        end;

        if (kickrank<0)
            kickrank=kickrank-1;
        end;

        if (max_dy<tol*als_tol_low)&&(kickrank>=0)
            kickrank=-1;
        end;

        if (max_dy<tol)&&(kickrank<=-als_iters)
            last_sweep=true;
            kickrank=0;
        end;

        if (order_index>numel(cur_order)) % New global sweep
            cur_order = block_order;
            order_index = 1;
            %residue
            if (last_sweep)
                cur_order = d-1;
            end;

            max_dy = 0;

%             dy_old = dy;
            swp = swp+1;
        end;

        dir = sign(cur_order(order_index));
        i = i+dir;
    end;
end

if (tailranks(1))
    cry{2} = reshape(cry{1}, n(1), ry(2))*reshape(cry{2}, ry(2), n(2)*ry(3));
    cry{2} = reshape(cry{2}, n(1), n(2), ry(3));    
end;
if (tailranks(2))
    cry{d-1} = reshape(cry{d-1}, ry(d-1)*n(d-1), ry(d))*reshape(cry{d}, ry(d), n(d));
    cry{d-1} = reshape(cry{d-1}, ry(d-1), n(d-1), n(d));
    cry = cry(1:d-1);
end;
if (tailranks(1))
    cry = cry(2:end);
end;
y = cell2core(y, cry);

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
