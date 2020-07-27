function [B]=tt_broadcast(A, n_extra)

% Simillar to Numpy function 'numpy.broadcast_to'.
% this function broadcast TT tensor A at extra dimensions with specified
% mode sizes.

% return a new TT tensor with B.n=[A.n; n_extra], and B.d=A.d +length(n_extra), 
% containing prod(n_extra) copies of A. 

% n_extra: vector of integers, gives the mode sizes of extra dimensions.
% length(n_extra): number of extra dimensions.

sn=size(n_extra);
if sn(2)>sn(1)
    n_extra=n_extra';
end
extra_dim= length(n_extra);

B=tt_ones(n_extra,extra_dim);
 
B.ps=[  A.ps(1:end) ; A.ps(end)+cumsum(n_extra)];
B.core= [A.core; B.core];

B.n = [A.n; B.n];
B.r = [A.r; B.r(2:end)];
B.d = A.d+B.d; 
end

