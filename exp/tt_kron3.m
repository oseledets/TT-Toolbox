
function [ttmm1] = tt_kron3(ttm1)
%[TTM11]=TT_KRON3(MAT,N2)
%Compute pseudo-Kronecker representation of the matrix
d = size(ttm1, 1);
ttmm1 = cell(d,1);

ns1 = zeros(d,1);
ns2 = zeros(d,1);

for q=1:d
    n1 = size(ttm1{q},1);
    m = size(ttm1{q},2);
   % n2 = size(ttm2{q},2);
    n2=m;
    r1 = size(ttm1{q}, 3);
    r2 = size(ttm1{q}, 4);
    ns1(q)=n1;
    ns2(q)=n2;
    
    I = eye(n2);
        ttmm1{q}=reshape(ttm1{q}, n1, m*r1*r2);
    ttmm1{q}=kron(I,ttmm1{q});
    ttmm1{q}=reshape(ttmm1{q}, n1, n2, m, r1, r2, n2);
    ttmm1{q}=permute(ttmm1{q}, [1 2 3 6 4 5]);
    ttmm1{q}=reshape(ttmm1{q}, n1*n2, m*n2, r1, r2);

end

