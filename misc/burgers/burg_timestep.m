function [unew]=burg_timestep(uold, f, tau, B,A,M)

n = size(uold,1);

I = speye(n);
rhs = f + uold/tau - B*uold*0.5 - (M*uold).*(A*uold)*0.25;

B_ts = I/tau + B*0.5 + spdiags(M*uold, 0, n, n)*A*0.25 + spdiags(A*uold, 0, n, n)*M*0.25;

% unew = newton_solve(B_ts,M*0.5,A*0.5,rhs,AMh*0.25,1e-10,uold);
unew = cg_solve(B_ts,M*0.5,A*0.5,rhs,1e-10,uold);

end