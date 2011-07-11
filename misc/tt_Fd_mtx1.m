% FUNCTION [M] = TT_FD_MTX1(TT_A, BOUND1, BOUND2)
% Generates finite-difference 1D matrix of \nabla(tt_a \nabla) in QTTM format.
% tt_a - in QTT format.
% bound1,2 - specify boundary conditions as x=0 and x=1, respectively:
%     0 - Dirichlet,
%     1 - Neumann
% eps - QTT compression tolerance

function [M] = tt_Fd_mtx1(tt_a, bound1, bound2, eps)

d = size(tt_a, 1);
M = tt_diag(tt_a);
tt_shf_a = tt_mv(tt_shf(d), tt_a);
tt_last_elem_a = tt_a;
for i=1:d
    tt_last_elem_a{i}(1,:,:)=0;
end;
M = ttm_add(M, tt_diag(tt_shf_a));
M = ttm_add(M, tt_diag(tt_last_elem_a));
M = tt_mat_compr(M, eps);
M_shf = tt_shf_diag(tt_a, 0);
M = ttm_add(M, tt_scal(M_shf, -1));
M = ttm_add(M, tt_scal(tt_transp(M_shf), -1));
M = tt_mat_compr(M, eps);
M = tt_scal(M, (2^d+1)^2);

end