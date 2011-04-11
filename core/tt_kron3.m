function [res] = tt_kron3(ttm1, ttm2)
%[res] = tt_kron3(ttm1, ttm2)
%Computes pseudo-Kronecker (i1j1,i2j2) product of the matrices.
%Notice that the "true" kronecker product is ttm2 \otimes ttm1:
%Here the indices of ttm1 vary faster than ttm2

d = size(ttm1, 1);
ns1 = tt_size(ttm1);
ms1 = tt_size(tt_transp(ttm1));
ns2 = tt_size(ttm2);
ms2 = tt_size(tt_transp(ttm2));
res = cell(d,1);
ttm1{1}=reshape(ttm1{1}, ns1(1), ms1(1), 1, size(ttm1{1},3));
ttm2{1}=reshape(ttm2{1}, ns2(1), ms2(1), 1, size(ttm2{1},3));

for q=1:d
    r11 = size(ttm1{q},3);
    r12 = size(ttm1{q},4);
    r21 = size(ttm2{q},3);
    r22 = size(ttm2{q},4);
    
    core1 = reshape(ttm1{q}, ns1(q), ms1(q)*r11*r12);
    core2 = reshape(ttm2{q}, ns2(q), ms2(q)*r21*r22);
    
    core = kron(core2, core1); % size ns1*ns2, ms1(q)*r11*r12*ms2(q)*r21*r22
    
    core = reshape(core, ns1(q),ns2(q), ms1(q),r11,r12,ms2(q),r21,r22);
    core = permute(core, [[1 2] [3 6] [4 7] [5 8]]);
    res{q}=reshape(core, ns1(q)*ns2(q), ms1(q)*ms2(q), r11*r21, r12*r22);
end;

res{1}=reshape(res{1}, size(res{1},1), size(res{1},2), size(res{1},4));

end