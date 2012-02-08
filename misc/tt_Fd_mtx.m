function [ttm] = tt_Fd_mtx(d_phys, tt_a, bound1, bound2, eps)
%Generates finite-difference diffusion matrix in QTT
%   [TTM] = TT_FD_MTX(D_PHYS, TT_A, BOUND1, BOUND2, PEPS) Generates 
%   finite-difference diffusion matrix in QTTM from the coefficient TT_A 
%   in QTT, D_PHYS - physical dimension of the problem, BOUND1,2 - 
%   boundary conditions at 0 and 1:
%       0 - Dirichlet,
%       1 - Neumann
%   EPS is the accuracy of the approximation. This function uses TT1.0
%   format.
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



d = size(tt_a, 1);
d1 = round(d/d_phys);
dims = ones(1,d_phys-1);
for i=1:d_phys-1
    dims(i) = d1*i;
end;
ranks = tt_ranks(tt_a);
phys_ranks = ranks(dims);
phys_j = ones(1, d_phys-1);
first_comp=1;
% ttm = tt_zeros(d,2);

while (phys_j(d_phys-1)<=phys_ranks(d_phys-1))
    
    for q=1:d_phys
        M = cell(d, 1);
        cur_a = cell(d1,1);

        for k=1:d1-1
            cur_a{k}=tt_a{k};
        end;
        cur_a{d1} = tt_a{d1}(:,:,phys_j(1));
        if (q==1)
            cur_a = tt_Fd_mtx1(cur_a, bound1, bound2,eps);
        else
            cur_a = tt_diag_av1(cur_a, eps);
        end;
        for k=1:d1-1
            M{k}=cur_a{k};
        end;
        M{d1}(1:2,1:2,1:size(cur_a{d1}, 3),1)=reshape(cur_a{d1}, 2, 2, size(cur_a{d1}, 3), 1);

        
        for p=2:d_phys-1
            cur_a{1} = permute(tt_a{d1*(p-1)+1}(:,phys_j(p-1),:), [1 3 2]);
            for k=2:d1-1
                cur_a{k}=tt_a{d1*(p-1)+k};
            end;
            cur_a{d1} = tt_a{d1*p}(:,:,phys_j(p));
            
            if (p==q)
                cur_a = tt_Fd_mtx1(cur_a, bound1, bound2,eps);
            else
                cur_a = tt_diag_av1(cur_a, eps);
            end;
            
            M{d1*(p-1)+1}(1:2,1:2,1,1:size(cur_a{1},3)) = reshape(cur_a{1}, 2, 2, 1, size(cur_a{1},3));
            for k=2:d1-1
                M{d1*(p-1)+k}(1:2,1:2,1:size(cur_a{k},3),1:size(cur_a{k},4))=cur_a{k};
            end;
            M{d1*p}(1:2,1:2,1:size(cur_a{d1},3), 1) = reshape(cur_a{d1}, 2, 2, size(cur_a{d1},3), 1);
        end;
        
        cur_a{1} = permute(tt_a{d1*(d_phys-1)+1}(:,phys_j(d_phys-1),:), [1 3 2]);
        for k=2:d1-1
            cur_a{k}=tt_a{d1*(d_phys-1)+k};
        end;
        cur_a{d1} = tt_a{d1*d_phys}(:,:,1);
        
        if (q==d_phys)
            cur_a = tt_Fd_mtx1(cur_a, bound1, bound2,eps);
        else
            cur_a = tt_diag_av1(cur_a, eps);
        end;
        
        M{d1*(d_phys-1)+1}(1:2,1:2,1,1:size(cur_a{1},3)) = reshape(cur_a{1}, 2, 2, 1, size(cur_a{1},3));
        for k=2:d1-1
            M{d1*(d_phys-1)+k}(1:2,1:2,1:size(cur_a{k},3),1:size(cur_a{k},4))=cur_a{k};
        end;
        M{d1*d_phys}(1:2,1:2,1:size(cur_a{d1},3), 1) = reshape(cur_a{d1}, 2, 2, size(cur_a{d1},3), 1);
        
        if (first_comp==1)
            ttm = M;
            first_comp=0;
        else
            ttm = ttm_add(ttm, M);
            ttm = tt_mat_compr(ttm, eps);
        end;
    end;
    
    phys_j(1)=phys_j(1)+1;
    for i=1:d_phys-2
        if (phys_j(i)>phys_ranks(i))
            phys_j(i)=1;
            phys_j(i+1)=phys_j(i+1)+1;
        end;
    end;
end;

end