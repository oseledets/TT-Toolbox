function [u, v] = cross2d_mat(a, eps)
%Classical maxvol-based cross
%       [U, V] = CROSS2D_MAT(A, EPS)
%       Computes the maxvol-based cross with requested accuracy EPS
%       A is a given matrix
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
eps0 = 1e-3; %Trunc parameter
[n, m] = size(a);
u = randn(n, r0);
v = randn(m, r0);

Phi = zeros(r0); 
er = 2 * eps;
ru = r0;
rv = r0;
ru_add = 0; 
rv_add = 0;
while ( er > eps )
    % Compute Phi
    [u, ru_mat] = qr(u, 0);
    [v, rv_mat] = qr(v, 0); 
    
    Phi = ru_mat(:, 1:ru) * Phi * (rv_mat(:, 1:rv)).'; 

    ru = ru + ru_add;
    rv = rv + rv_add;
    
    indu = maxvol2(u);
    indv = maxvol2(v); 

    sbm = a(indu, indv);
    sbmu = u(indu, :);
    sbmv = v(indv, :);
    Phi_new = sbmu \ sbm / (sbmv.');
    %Compute the error
    er = norm(Phi_new - Phi,'fro')/norm(Phi_new, 'fro');
    fprintf('ru = %d, rv = %d, er = %3.1e \n', ru, rv, er);
    Phi = Phi_new;
    if ( er > eps )
        %Compute new columns
        
        uadd = a(:, indv);
        vadd = a(indu, :).';
        
        %Compute truncated addition
        [uadd, sadd, dmp] = svd(uadd, 'econ');
        sadd = diag(sadd);
        ru_add = my_chop2(sadd, norm(sadd) * eps0);
        ru_add = min(ru + ru_add, n) - ru;
        uadd = uadd(:, 1:ru_add);
        [vadd, sadd, dmp] = svd(vadd, 'econ');
        sadd = diag(sadd);
        rv_add = my_chop2(sadd, norm(sadd) * eps0);
        rv_add = min(rv + rv_add, m) - rv;
        vadd = vadd(:, 1:rv_add);
        u = [u, uadd];
        v = [v, vadd];
    end
end
[u0, s0, v0] = svd(Phi, 'econ');
s0 = diag(s0); r0 = my_chop2(s0, eps * norm(s0));
u0 = u0(:, 1:r0); s0 = s0(1:r0); v0 = v0(:, 1:r0);
u = u * u0;
v = v * v0;
u = u * diag(s0);
fprintf('Done with rank = %d\n', r0);
end