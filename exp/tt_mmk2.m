function [ttm12] = tt_mmk2(ttm1, ttm2, eps, max_r, max_swp)
%DMRG/Krylov matrix-by-matrix multiplication 
%   [TTM12] = TT_MMK2(TTM1, TTM2, EPS, [MAX_R, MAX_SWP]) 
%   Computes DMRG/Krylov matrix-by-matrix multiplication with accuracy EPS.
%   It just calls the MVK2 code. Uses TT1.0, so it syntax will change in
%   the next release
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

d = size(ttm1, 1);

exists_max_r=1;
if ((nargin<4)||(isempty(max_r)))
    exists_max_r=0;
end;
exists_max_swp=1;
if ((nargin<5)||(isempty(max_swp)))
    exists_max_swp=0;
end;

ttmm1 = cell(d,1);

ns1 = zeros(d,1);
ns2 = zeros(d,1);

for q=1:d
    n1 = size(ttm1{q},1);
    m = size(ttm1{q},2);
    n2 = size(ttm2{q},2);
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
%     ttmm1{q} = zeros(n1*n2, m*n2, size(ttm1{q}, 3), size(ttm1{q}, 4));
    
%     ttmm1{q}(1,1,:,:)=ttm1{q}(1,1,:,:);
%     ttmm1{q}(2,1,:,:)=ttm1{q}(1,2,:,:);
%     ttmm1{q}(1,2,:,:)=ttm1{q}(2,1,:,:);
%     ttmm1{q}(2,2,:,:)=ttm1{q}(2,2,:,:);
%     
%     ttmm1{q}(3,3,:,:)=ttm1{q}(1,1,:,:);
%     ttmm1{q}(4,3,:,:)=ttm1{q}(1,2,:,:);
%     ttmm1{q}(3,4,:,:)=ttm1{q}(2,1,:,:);
%     ttmm1{q}(4,4,:,:)=ttm1{q}(2,2,:,:);    
end;

ttmm2 = tt_mat_to_vec(ttm2);

if (exists_max_r)
    if (exists_max_swp)
        ttm12 = tt_mvk2(ttmm1, ttmm2, eps, max_r, max_swp);
    else
        ttm12 = tt_mvk2(ttmm1, ttmm2, eps, max_r);
    end;
else
    if (exists_max_swp)
        ttm12 = tt_mvk2(ttmm1, ttmm2, eps, [], max_swp);
    else
        ttm12 = tt_mvk2(ttmm1, ttmm2, eps);
    end;
end;

ttm12 = tt_vec_to_mat(ttm12, ns1, ns2);

end