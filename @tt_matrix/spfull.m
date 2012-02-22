function [B]=spfull(A, tol)

if (nargin<2)||(isempty(tol))
    tol = 1e-14;
end;

tA = A.tt;
d = tA.d;
r = tA.r;
n = A.n;
m = A.m;

B0 = cell(1,1);
B0{1} = 1;

for i=1:d
    curB = tA{i};
    curB = reshape(curB, r(i), n(i), m(i), r(i+1));
    B = cell(r(i+1),1);
    for k2=1:r(i+1) % Right index - hang it for now
        B{k2} = sparse(prod(n(1:i)), prod(m(1:i)));
        for k1=1:r(i) % Left index - convolve with the current part of B
            curM = curB(k1, :, :, k2);
            curM = reshape(curM, n(i), m(i));
            maxM = max(abs(curM(:)));
            curM(abs(curM)<tol*maxM)=0;
            curM = sparse(curM);
            B{k2} = B{k2} + kron(curM, B0{k1});
        end;
    end;
    B0 = B;
end;

B = B{1};
end