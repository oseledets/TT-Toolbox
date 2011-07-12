function [z]=conv3full(a, x, y)
% function [z]=conv3full(a, x, y)
% Convolves the 3D tt_array with two full vectors x,y over
% 2nd and 3d modes (use permute2(a,order) if necessary), respectively
% (quadratic convolution).

nsa = a.ns;
d = size(nsa,1);
dp = size(nsa,2);
if (dp~=3)
    error('High-order convolutions are not implemented yet. Only Q_{ijk} x_j y_k now.\n');
end;

tta = a.tt;
ra = tta.r;
psa = tta.ps;
cra = tta.core;

% First mode is not trivial.
% Convolve both x and y
N2 = prod(nsa(2:d,2));
N3 = prod(nsa(2:d,3));

z = reshape(x, nsa(1,2), N2);
a1 = cra(psa(1):psa(2)-1);
a1 = reshape(a1, [ra(1), nsa(1,1), nsa(1,2), nsa(1,3), ra(2)]);
a1 = permute(a1, [2,5,4,1,3]);
a1 = reshape(a1, nsa(1,1)*ra(2)*nsa(1,3), ra(1)*nsa(1,2));
z = a1*z; % size n1*ra2*n3, N2
z = reshape(z, nsa(1,1), ra(2), nsa(1,3), N2);
z = permute(z, [1,4,2,3]);
z = reshape(z, nsa(1,1)*N2*ra(2), nsa(1,3));

cur_y = reshape(y, nsa(1,3), N3);
z = z*cur_y; % size n1*N2*ra2, N3
N1 = nsa(1,1);
N2 = round(N2/nsa(2,2));
N3 = round(N3/nsa(2,3));
z = reshape(z, N1, nsa(2,2), N2, ra(2), nsa(2,3), N3);
z = permute(z, [4 2 5  1 3 6]);
z = reshape(z, ra(2)*nsa(2,2)*nsa(2,3), N1*N2*N3);

% Successively convolve middle cores
for i=2:d-1
    a1 = cra(psa(i):psa(i+1)-1);
    a1 = reshape(a1, [ra(i), nsa(i,1), nsa(i,2), nsa(i,3), ra(i+1)]);    
    a1 = permute(a1, [2,5,1,3,4]);
    a1 = reshape(a1, nsa(i,1)*ra(i+1), ra(i)*nsa(i,2)*nsa(i,3));
    
    z = a1*z; % size n1*ra2, N1*N2*N3
    N2 = round(N2/nsa(i+1,2));
    N3 = round(N3/nsa(i+1,3));
    z = reshape(z, nsa(i,1), ra(i+1), N1, nsa(i+1,2), N2, nsa(i+1,3), N3);
    z = permute(z, [2 4 6  3 1 5 7]);
    N1 = N1*nsa(i,1);
    z = reshape(z, ra(i+1)*nsa(i+1,2)*nsa(i+1,3), N1*N2*N3);
end;

% The last mode
% z = z[r(d-1)*nsa(d,2)*nsa(d,3), prod(nsa(1:d-1,1)]
% We are to convolve the first dimension
a1 = cra(psa(d):psa(d+1)-1);
a1 = reshape(a1, [ra(d), nsa(d,1), nsa(d,2), nsa(d,3), ra(d+1)]);
a1 = permute(a1, [2,5,1,3,4]);
a1 = reshape(a1, nsa(d,1)*ra(d+1), ra(d)*nsa(d,2)*nsa(d,3));

z = a1*z; % size nsa(d,1),N1;
z = reshape(z.', N1*nsa(d,1));
end
