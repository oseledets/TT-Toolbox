function test_tt_gmres_complex()
% From the toolbox root: setup; addpath('tests'); test_tt_gmres_complex
% Four distinct complex eigenvalues exercise complex Arnoldi projections.
lambda = [1+1i; 2-1i; 3+2i; 4-0.5i];
exact = [1; -1i; 2+1i; 1-2i];
b = lambda.*exact;
A = diag(tt_tensor(reshape(lambda,[2 2]),1e-14));
rhs = tt_tensor(reshape(b,[2 2]),1e-14);
[x,~] = tt_gmres(core(A),core(rhs),1e-10,1,4,1e-14,1e-14,[],[],[],[],0);
field = reshape(x{1}*x{2}.',[],1);
assert(norm(lambda.*field-b)/norm(b) < 1e-9);
assert(norm(field-exact)/norm(exact) < 1e-9);
end
