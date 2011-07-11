function [AMh] = newton_prepare(A, M)
% Prepares M_{\tilde{i}, i} \otimes A^T(\hat{i}, i) for the Hessian part;
% Involves BAD MATLAB cycles

n = size(A,1);

AMh = zeros(n,n^2);
pre_AMh = kron(M, A); % indices i*i, \tilde{i}*\hat{i}

for i=1:n
    AMh(i,:) = full(pre_AMh(i+(i-1)*n, :));
end;

AMh = AMh';

end