function [unew]=burg_timestep(uold, f, tau, B,A,M, AMh)

n = size(uold,1);

I = speye(n);
rhs = f + uold/tau - B*uold - (M*uold).*(A*uold);

B_ts = I/tau + B + spdiags(M*uold, 0, n, n)*A + spdiags(A*uold, 0, n, n)*M;

unew = newton_solve(B_ts,M,A,rhs,AMh,1e-10,uold);

end