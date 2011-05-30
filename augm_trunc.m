function [X,Y]=augm_trunc(A, u, r, beta)

n1 = size(u,1);
n2 = size(u,2);
maxit=50;

% X = randn(n1,r);
% Y = randn(n2,r);
[X,S,Y]=svd(u,'econ');
X = X(:,1:r);
Y = Y(:,1:r)*S(1:r,1:r);
Y=conj(Y);

mu = ones(n1*n2,1);
% rhs = tt_mat_full_vec(A,reshape(u, n1*n2, 1));
rhs = A*reshape(u, n1*n2, 1);
for i=1:maxit
%     C = A'*A + spdiags(mu, [0], n1*n2, n1*n2);
%     C = tt_mm(tt_transp(A), A);
%     C{1}(:,:,size(C{1},3)+1)=beta*eye(n1);
%     C{2}(:,:,size(C{2},3)+1)=eye(n2);
%     C = tt_mat_compr(C, 1e-12);
%     iC = tt_minres_selfprec(C, 1e-1, 1e-3, 15, 'right');    
    C = A'*A + beta*speye(size(A,2));
%     rhs_new = tt_mat_full_vec(tt_transp(A),rhs) + beta*reshape(X*Y.', n1*n2, 1) + mu;
    rhs_new = A'*rhs + beta*reshape(X*Y.', n1*n2, 1) + mu;
%     cur_u = gmres(@(v)tt_mat_full_vec(C,v), rhs_new, 50, 1e-10, 50, @(v2)tt_mat_full_vec(iC,v2));
    cur_u = C \ rhs_new;
    
%     Mx = Y'*Y;
%     rhs_x = Y'*reshape(-mu+cur_u,n1,n2)';
%     X = pcg(kron(Mx, eye(n1)), reshape(rhs_x, r*n1, 1), 1e-8, 1000);
% %     X = Mx \ rhs_x;
%     X = reshape(X, r, n1)';
%     
%     My = X'*X;
%     rhs_y = X'*reshape(-mu+cur_u,n1,n2);
% %     Y = My \ rhs_y;
%     Y = pcg(kron(My, eye(n2)), reshape(rhs_y, r*n2, 1), 1e-8, 1000);
%     Y = reshape(Y, r, n2)';    
%     Y = Y';    
    cur_u = reshape(cur_u, n1, n2);    
    [X,S,Y]=svd(cur_u,'econ');
    X = X(:,1:r);
    Y = Y(:,1:r)*S(1:r,1:r);
    Y=conj(Y);
    
    err = reshape(cur_u, n1, n2) - X*Y.';
    mu = mu - beta*reshape(err, n1*n2, 1);
%     fprintf('iter: %d, er0=%3.3e, mu=%3.3e, resid = %3.3e\n', i, norm(err(:))/norm(cur_u(:)), norm(mu), norm(tt_mat_full_vec(A,reshape(X*Y.', n1*n2, 1))-rhs)/norm(rhs));
    fprintf('iter: %d, er0=%3.3e, mu=%3.3e, resid = %3.3e\n', i, norm(err(:))/norm(cur_u(:)), norm(mu), norm(A*reshape(X*Y.', n1*n2, 1)-rhs)/norm(rhs));
end;

end