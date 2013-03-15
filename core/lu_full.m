function [ind] = lu_full(a)
%LU with full pivoting
%    [U,V,IND] = LU_FULL(A)
b = a;
[n, r] = size(b);
ind = zeros(r, 1);
for i = 1:r
    [mx0,big_ind] = max(abs(b(:)));
    [i0,j0] = ind2sub([n,r],big_ind);
    b = b - b(:, j0) * b(i0, :) / b(i0,j0);
    ind(i) = i0;
end

