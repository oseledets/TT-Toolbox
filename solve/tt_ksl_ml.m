function [x]=tt_ksl_ml(x0, A, y, tau)
% function [x]=tt_ksl_ml(x0, A, y, tau)
% Linear KSL time scheme in the TT format
% Integrates the system dx/dt + Ax = y on t=[0, tau].
% Uses the Crank-Nicolson as the local solver:
% the matrix A is supposed to be time-independent, and R.H.S. y is given at
% the middle time point.

n = x0.n;
d = x0.d;
rx = x0.r;
ra = A.tt;
ra = ra.r;
ry = y.r;
% ry2 = y2.r;

crx = core2cell(x0);
cry = core2cell(y);
% cry2 = core2cell(y2);
cra = core2cell(A);

phia = cell(d+1,1); phia{1}=1; phia{d+1}=1;
phiy = cell(d+1,1); phiy{1}=1; phiy{d+1}=1;
% phiy2 = cell(d+1,1); phiy2{1}=1; phiy2{d+1}=1;

% tau = tau/2;

% 
%%%% We do not need actually the backward sweep. It only makes the scheme
%%%% unstable.
%

% % orth, phi
% for i=1:d-1
%     cr = reshape(crx{i}, rx(i)*n(i), rx(i+1));
%     [cr,rv]=qr(cr, 0);
%     cr2 = reshape(crx{i+1}, rx(i+1), n(i+1)*rx(i+2));
%     cr2 = rv*cr2;
%     rx(i+1) = size(cr, 2);
%     crx{i} = reshape(cr, rx(i), n(i), rx(i+1));
%     crx{i+1} = reshape(cr2, rx(i+1), n(i+1), rx(i+2));
%     
%     phia{i+1} = compute_next_Phi(phia{i}, crx{i}, cra{i}, crx{i}, 'lr');
%     phiy{i+1} = compute_next_Phi(phiy{i}, crx{i}, [], cry{i}, 'lr');    
%     phiy2{i+1} = compute_next_Phi(phiy2{i}, crx{i}, [], cry2{i}, 'lr');    
% end;


% orth, phi
for i=d:-1:2
    cr = reshape(crx{i}, rx(i), n(i)*rx(i+1));
    [cr,rv]=qr(cr.', 0);
    cr2 = reshape(crx{i-1}, rx(i-1)*n(i-1), rx(i));
    cr2 = cr2*rv.';
    rx(i) = size(cr, 2);
    crx{i} = reshape(cr.', rx(i), n(i), rx(i+1));
    crx{i-1} = reshape(cr2, rx(i-1), n(i-1), rx(i));
    
    phia{i} = compute_next_Phi(phia{i+1}, crx{i}, cra{i}, crx{i}, 'rl');
    phiy{i} = compute_next_Phi(phiy{i+1}, crx{i}, [], cry{i}, 'rl');    
%     phiy2{i} = compute_next_Phi(phiy2{i+1}, crx{i}, [], cry2{i}, 'rl');    
end;

% forward: propagation1
for i=1:d
    %K
    
    Phi1 = phia{i}; Phi2 = phia{i+1};
    % Phi1: rx'1, rx1, ra1, or rx'1, ry1
    % Phi2: rx2, ra2, rx'2, or ry'2, rx2
    A1 = cra{i}; y1 = cry{i};
    % RHS - rewrite it in accordance with new index ordering
    rhs = phiy{i}; % rx'1, ry1
    y1 = reshape(y1, ry(i), n(i)*ry(i+1));
    rhs = rhs*y1;
    rhs = reshape(rhs, rx(i)*n(i), ry(i+1));
    rhs = rhs*phiy{i+1}; 
    rhs = reshape(rhs, rx(i)*n(i)*rx(i+1),1);
%     rhs2 = phiy2{i}; % rx'1, ry1
%     y1 = reshape(cry2{i}, ry2(i), n(i)*ry2(i+1));
%     rhs2 = rhs2*y1;
%     rhs2 = reshape(rhs2, rx(i)*n(i), ry2(i+1));
%     rhs2 = rhs2*phiy2{i+1}; 
%     rhs2 = reshape(rhs2, rx(i)*n(i)*rx(i+1),1);
%     
    % sol_prev
    sol_prev = reshape(crx{i}, rx(i)*n(i)*rx(i+1), 1);

%     if (rx(i)*n(i)*rx(i+1)<max_full_size) % Full solution
        %      |     |    |
        % B = Phi1 - A1 - Phi2
        %      |     |    |
        B = reshape(Phi1, rx(i)*rx(i), ra(i));
        B = B*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
        B = reshape(B, rx(i), rx(i), n(i), n(i)*ra(i+1));
        B = permute(B, [1, 3, 2, 4]);
        B = reshape(B, rx(i)*n(i)*rx(i)*n(i), ra(i+1));
        B = B*reshape(permute(Phi2, [2, 3, 1]), ra(i+1), rx(i+1)*rx(i+1));
        B = reshape(B, rx(i)*n(i), rx(i)*n(i), rx(i+1), rx(i+1));
        B = permute(B, [1, 3, 2, 4]);
        B = reshape(B, rx(i)*n(i)*rx(i+1), rx(i)*n(i)*rx(i+1));
        
    CNm = eye(rx(i)*n(i)*rx(i+1)) - tau*0.5*B;
    CNp = eye(rx(i)*n(i)*rx(i+1)) + tau*0.5*B;

%     sol1 = CNm*sol_prev + tau*(rhs+rhs2)*0.5;
    sol1 = CNm*sol_prev + tau*rhs;
    sol = CNp \ sol1;

    if (i<d) % we have S and V(implicitly)
        
        %S
        sol = reshape(sol, rx(i)*n(i), rx(i+1));
        [U,S_prev,V_prev] = svd(sol, 'econ');
        S_prev = S_prev*V_prev';
        rnew = size(U,2);
        % Project matrix, rhs.
        Bs = reshape(B, rx(i)*n(i), rx(i+1)*rx(i)*n(i)*rx(i+1));
        Bs = U'*Bs;
        Bs = reshape(Bs, rnew*rx(i+1), rx(i)*n(i), rx(i+1));
        Bs = permute(Bs, [1,3,2]);
        Bs = reshape(Bs, rnew*rx(i+1)*rx(i+1), rx(i)*n(i));
        Bs = Bs*U;
        Bs = reshape(Bs, rnew*rx(i+1), rx(i+1), rnew);
        Bs = permute(Bs, [1,3,2]);
        Bs = reshape(Bs, rnew*rx(i+1), rnew*rx(i+1));
        
        rhss = reshape(rhs, rx(i)*n(i), rx(i+1));
        rhss = U'*rhss;
        rhss = reshape(rhss, rnew*rx(i+1), 1);
%         rhss2 = reshape(rhs2, rx(i)*n(i), rx(i+1));
%         rhss2 = U'*rhss2;
%         rhss2 = reshape(rhss2, rnew*rx(i+1), 1);
        
        CNm = eye(rnew*rx(i+1)) + tau*0.5*Bs; % !!!! backward step
        CNp = eye(rnew*rx(i+1)) - tau*0.5*Bs; % !!!! backward step
        
        S_prev = reshape(S_prev, rnew*rx(i+1), 1);
        S1 = CNm*S_prev - tau*rhss;
%         S1 = CNm*S_prev - tau*(rhss+rhss2)*0.5;        
        S = CNp \ S1;
        
        % Stuff back and recompute phis
        S = reshape(S, rnew, rx(i+1));
        cr2 = reshape(crx{i+1}, rx(i+1), n(i+1)*rx(i+2));
        cr2 = S*cr2;
        rx(i+1) = rnew;
        crx{i+1} = reshape(cr2, rx(i+1), n(i+1), rx(i+2));
        crx{i} = reshape(U, rx(i), n(i), rx(i+1));
        
        phia{i+1} = compute_next_Phi(phia{i}, crx{i}, cra{i}, crx{i}, 'lr');
        phiy{i+1} = compute_next_Phi(phiy{i}, crx{i}, [], cry{i}, 'lr');
%         phiy2{i+1} = compute_next_Phi(phiy2{i}, crx{i}, [], cry2{i}, 'lr');
        
        % V is done implicitly via recursion
    else
        % K is V, stuff back
        crx{i} = reshape(sol, rx(i), n(i), rx(i+1));
    end;
end;

%%% vvv  Is not needed  vvv
%%%
% % backward: propagation2
% for i=d:-1:1
%     %V
%     
%     Phi1 = phia{i}; Phi2 = phia{i+1};
%     % Phi1: rx'1, rx1, ra1, or rx'1, ry1
%     % Phi2: rx2, ra2, rx'2, or ry'2, rx2
%     A1 = cra{i}; y1 = cry{i};
%     % RHS - rewrite it in accordance with new index ordering
%     rhs = phiy{i}; % rx'1, ry1
%     y1 = reshape(y1, ry(i), n(i)*ry(i+1));
%     rhs = rhs*y1;
%     rhs = reshape(rhs, rx(i)*n(i), ry(i+1));
%     rhs = rhs*phiy{i+1}; 
%     rhs = reshape(rhs, rx(i)*n(i)*rx(i+1),1);
%     rhs2 = phiy2{i}; % rx'1, ry1
%     y1 = reshape(cry2{i}, ry2(i), n(i)*ry2(i+1));
%     rhs2 = rhs2*y1;
%     rhs2 = reshape(rhs2, rx(i)*n(i), ry2(i+1));
%     rhs2 = rhs2*phiy2{i+1}; 
%     rhs2 = reshape(rhs2, rx(i)*n(i)*rx(i+1),1);
%     
%     % sol_prev
%     sol_prev = reshape(crx{i}, rx(i)*n(i)*rx(i+1), 1);
% 
% %     if (rx(i)*n(i)*rx(i+1)<max_full_size) % Full solution
%         %      |     |    |
%         % B = Phi1 - A1 - Phi2
%         %      |     |    |
%         B = reshape(Phi1, rx(i)*rx(i), ra(i));
%         B = B*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
%         B = reshape(B, rx(i), rx(i), n(i), n(i)*ra(i+1));
%         B = permute(B, [1, 3, 2, 4]);
%         B = reshape(B, rx(i)*n(i)*rx(i)*n(i), ra(i+1));
%         B = B*reshape(permute(Phi2, [2, 3, 1]), ra(i+1), rx(i+1)*rx(i+1));
%         B = reshape(B, rx(i)*n(i), rx(i)*n(i), rx(i+1), rx(i+1));
%         B = permute(B, [1, 3, 2, 4]);
%         B = reshape(B, rx(i)*n(i)*rx(i+1), rx(i)*n(i)*rx(i+1));
%         
%     CNm = eye(rx(i)*n(i)*rx(i+1)) - tau*0.5*B;
%     CNp = eye(rx(i)*n(i)*rx(i+1)) + tau*0.5*B;
% 
%     sol1 = CNm*sol_prev + tau*(rhs+rhs2)*0.5;
%     sol = CNp \ sol1;
% 
%     if (i>1) % we have S and K(implicitly)
%         
%         %S
%         sol = reshape(sol, rx(i), n(i)*rx(i+1));
%         [V,S_prev,V_prev] = svd(sol.', 'econ');
%         S_prev = (S_prev*V_prev').';
%         rnew = size(V,2);
%         % Project matrix, rhs.
%         Bs = reshape(B, rx(i)*n(i)*rx(i+1)*rx(i), n(i)*rx(i+1));
%         Bs = Bs*V;
%         Bs = reshape(Bs, rx(i), n(i)*rx(i+1), rx(i)*rnew);
%         Bs = permute(Bs, [2,1,3]);
%         Bs = reshape(Bs, n(i)*rx(i+1), rx(i)*rx(i)*rnew);
%         Bs = V'*Bs;
%         Bs = reshape(Bs, rnew, rx(i), rx(i)*rnew);
%         Bs = permute(Bs, [2,1,3]);
%         Bs = reshape(Bs, rx(i)*rnew, rx(i)*rnew);
%         
%         rhss = reshape(rhs, rx(i), n(i)*rx(i+1));
%         rhss = rhss*conj(V);
%         rhss = reshape(rhss, rx(i)*rnew, 1);
%         rhss2 = reshape(rhs2, rx(i), n(i)*rx(i+1));
%         rhss2 = rhss2*conj(V);
%         rhss2 = reshape(rhss2, rx(i)*rnew, 1);
%         
%         CNm = eye(rx(i)*rnew) + tau*0.5*Bs; % !!!! backward step
%         CNp = eye(rx(i)*rnew) - tau*0.5*Bs; % !!!! backward step
%         
%         S_prev = reshape(S_prev, rx(i)*rnew, 1);
%         S1 = CNm*S_prev - tau*(rhss+rhss2)*0.5;
%         S = CNp \ S1;
%         
%         % Stuff back and recompute phis
%         S = reshape(S, rx(i), rnew);
%         cr2 = reshape(crx{i-1}, rx(i-1)*n(i-1), rx(i));
%         cr2 = cr2*S;
%         rx(i) = rnew;
%         crx{i-1} = reshape(cr2, rx(i-1), n(i-1), rx(i));
%         crx{i} = reshape(V.', rx(i), n(i), rx(i+1));
%         
%         phia{i} = compute_next_Phi(phia{i+1}, crx{i}, cra{i}, crx{i}, 'rl');
%         phiy{i} = compute_next_Phi(phiy{i+1}, crx{i}, [], cry{i}, 'rl');
%         phiy2{i} = compute_next_Phi(phiy2{i+1}, crx{i}, [], cry2{i}, 'rl');
%         
%         % K is done implicitly via recursion
%     else
%         % K is V, stuff back
%         crx{i} = reshape(sol, rx(i), n(i), rx(i+1));
%     end;
% end;

x = cell2core(tt_tensor, crx);
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
