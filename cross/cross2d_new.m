function [u, v] = cross2d_new(f, n, m, eps, varargin)
%Classical maxvol-based cross
%       [U, V] = CROSS2D(F, N, M, EPS)
%       Computes the maxvol-based cross with requested accuracy EPS
%       F is a function that computes the prescribed element, F(I, J) 
%       N, M are the sizes of the matrix
%
% TT-Toolbox 2.2, 2009-2013
% 
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
 
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
r0 = 2;
full_check = false;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'full_check'
            full_check = varargin{i+1};
        case 'r0'
            r0 = varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end


u = randn(n, r0);
v = randn(m, r0);

er = 2 * eps;



[u, dmp] = qr(u, 0);
[v, dmp] = qr(v, 0);
indu = maxvol2(u);
indv = maxvol2(v);
uu = u / u(indu, :);
vv = v / v(indv, :);
%u = u \ sbmu;
%v = v \ sbmv;
indu_add = indu;
indv_add = indv; 
if ( full_check )
    fmat = zeros(n, m);
    for i = 1:n
        for j = 1:m
            fmat(i, j) = f(i, j);
        end
    end
end
ru = numel(indu);
rv = numel(indv);
Phi = zeros(ru, rv);

for s1 = 1:ru
    for s2 = 1:rv
        Phi(s1, s2) = f(indu(s1), indv(s2));
    end
end
while ( er > eps )
    %Compute uadd, vadd
    rv_add = numel(indv_add);
    ru_add = numel(indu_add);
    uadd = zeros(n, rv_add);
    for i = 1:n
        for s = 1:rv_add
            uadd(i, s) = f(i, indv_add(s));
        end
    end
    uadd_appr = uu * Phi * vv(indv_add, :).';
    vadd = zeros(m, ru_add);
    for j = 1:m
        for s = 1:ru_add
            vadd(j, s) = f(indu_add(s), j);
        end
    end
    vadd_appr = uu(indu_add, :) * Phi * vv.';
    vadd_appr = vadd_appr.';
    er1 = norm(uadd_appr - uadd, 'fro')/norm(uadd, 'fro');
    er2 = norm(vadd_appr - vadd, 'fro')/norm(vadd, 'fro');
    er = max(er1, er2);
    %Compute the Schur complement 
    %uadd - u * (uadd()
    %[A11 A12] [u] = [f]
    %[A21 A22] [v] = [g]
    % u = A11^{-1}(f - A12 * v)
    % [A22 - (A21 * A11^{-1}) A12)] v = 
    comp_u = uadd - uu * uadd(indu, :);
    [comp_u, dmp] = qr(comp_u, 0);
    %[U Uadd] * PhiNew * [V Vadd]'
    %I think the current scheme gives  you a 
    %very simple idea to get the (next) pivots
    
    indu_add = maxvol2(comp_u);
    %comp_u = comp_u / comp_u(indu_add, :); 
    %
    %Our goal is to update UU^{-1}
    %S * wU = U
    %[u1 u2] * [A11 A12] = [v1 v2]
    %          [A21 A22]
    %u1 * A11 + u2 * A21 = v1
    %u1 * A12 + u2 * A22 = v2
    %u1 = v1 * A11^{-1}  - u2 * A21 * A11^{-1}
    %u2 * (A22 - A21 * A11^{-1} * A12) = v2 
    %v1 * A11^{-1} * A12 - u2 * A21 * A11^{-1} * A12 + u2 * A22 = v2
    %u2 * S = v2 - v1 * A11^{-1} * A12
    %u2 = v2 * S^{-1} - v1 * A11^{-1} * A12 * S^{-1} 
    %u1 = v1 * A11^{-1} - S' * S^{-1} * A21 * A11{-1} 
    %}
    %[A21 A22]    [v] 
    u2 = comp_u / comp_u(indu_add, :);
    u1 = uu - u2 * uu(indu_add, :);
    uu = [u1, u2];
    fu = max(abs(uu(:)));
    indu1 = [indu, indu_add];

    comp_v = vadd - vv * vadd(indv, :);
    [comp_v, dmp] = qr(comp_v, 0);
    indv_add = maxvol2(comp_v);

    v2 = comp_v / comp_v(indv_add, :);
    v1 = vv - v2 * vv(indv_add, :);
    vv = [v1, v2];
    fv = max(abs(vv(:)));
    er = er * max(fu, fv);
    indv1 = [indv, indv_add];
    ru = numel(indu1);
    rv = numel(indv1);
    Phi = zeros(ru, rv);
    for s1 = 1:ru
        for s2 = 1:rv
            Phi(s1, s2) = f(indu1(s1), indv1(s2)); 
        end
    end
    indu = indu1;
    indv = indv1;
    
    if ( full_check)
        appr = uu * Phi * vv.';
        fprintf('tr: %3.1e est: %3.1e fu: %3.1f fv: %3.1f r: %d \n', norm(appr - fmat,'fro') / norm(appr, 'fro'), er2, fu, fv, size(uu, 2));
    end
end
 u = uu * Phi;
 v = vv;
end
