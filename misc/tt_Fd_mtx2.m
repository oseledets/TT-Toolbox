function [ttm] = tt_Fd_mtx2(tt_a, bound1, bound2, eps)
%TT-representation of the diffusion matrix
%   [TTM] = TT_FD_MTX2(TT_A, BOUND1, BOUND2, EPS) Computes TT
%   representation of a simplest discretization of the diffusion operator
%   with operator given in the QTT-format (TT_A). 
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
n = zeros(1,d);
for q=1:d
    n(q)=size(tt_a{q}, 1);
end;

if (max(size(bound1))==1) bound1 = bound1*ones(1,d); end;
if (max(size(bound2))==1) bound2 = bound2*ones(1,d); end;

ranks_a = tt_ranks(tt_a);

for q=1:d
    cur_tt = cell(d,1);
    
    cur_tt{1}=zeros(n(1),n(1),ranks_a(1));
    for k=1:ranks_a(1)
            lp1 = diag(-1*ones(n(1)-1,1), [-1])+diag(2*ones(n(1), 1), [0])+diag(-1*ones(n(1)-1,1), [1]);
            Mp1 = diag((1/6)*ones(n(1)-1,1), [-1])+diag((4/6)*ones(n(1), 1), [0])+diag((1/6)*ones(n(1)-1,1), [1]);                        
            if (bound1(1)==1) lp1(1,1)=1; end;
            if (bound2(1)==1) lp1(n(1),n(1))=1; end;
            h = 1/(n(1)+1);
            if ((bound1(1)==1)||(bound2(1)==1)) h = 1/n(1); end;
            if ((bound1(1)==1)&&(bound2(1)==1)) h = 1/(n(1)-1); end;
            Mlp1 = diag(tt_a{1}(2:n(1),k), [-1])+diag(tt_a{1}(2:n(1),k), [1]);
            Mlp1 = Mlp1 + diag(tt_a{1}(1:n(1),k)+[tt_a{1}(2:n(1),k)' tt_a{1}(n(1),k)]', [0])*0.5;
        if (q==1)            
            cur_tt{1}(:,:,k) = Mlp1.*(lp1/h^2);
        else            
            cur_tt{1}(:,:,k) = Mlp1.*Mp1;
%             cur_tt{1}(:,:,k) = diag(tt_a{1}(1:n(1),k)+[tt_a{1}(2:n(1),k)' tt_a{1}(n(1),k)]', [0])*0.5;
        end;
    end;
    
    for p=2:d-1
        cur_tt{p}=zeros(n(p),n(p),ranks_a(p-1),ranks_a(p));
        for k1=1:ranks_a(p-1)
            for k2=1:ranks_a(p)
                    lp1 = diag(-1*ones(n(p)-1,1), [-1])+diag(2*ones(n(p), 1), [0])+diag(-1*ones(n(p)-1,1), [1]);
                    Mp1 = diag((1/6)*ones(n(1)-1,1), [-1])+diag((4/6)*ones(n(1), 1), [0])+diag((1/6)*ones(n(1)-1,1), [1]);                                        
                    if (bound1(p)==1) lp1(1,1)=1; end;
                    if (bound2(p)==1) lp1(n(p),n(p))=1; end;
                    h = 1/(n(p)+1);
                    if ((bound1(p)==1)||(bound2(p)==1)) h = 1/n(p); end;
                    if ((bound1(p)==1)&&(bound2(p)==1)) h = 1/(n(p)-1); end;
                    Mlp1 = diag(tt_a{p}(2:n(p),k1,k2), [-1])+diag(tt_a{p}(2:n(p),k1,k2), [1]);
                    Mlp1 = Mlp1 + diag(tt_a{p}(1:n(p),k1,k2)+[tt_a{p}(2:n(p),k1,k2)' tt_a{p}(n(p),k1,k2)]', [0])*0.5;
                if (p==q)                           
                    cur_tt{p}(:,:,k1,k2) = Mlp1.*(lp1/h^2);             
                else
                    cur_tt{p}(:,:,k1,k2) = Mlp1.*Mp1;
%                     cur_tt{p}(:,:,k1,k2) = diag(tt_a{p}(1:n(p),k1,k2)+[tt_a{p}(2:n(p),k1,k2)' tt_a{p}(n(p),k1,k2)]', [0])*0.5;
                end;
            end;
        end;
    end;
    
    cur_tt{d}=zeros(n(d),n(d),ranks_a(d-1));
    for k=1:ranks_a(d-1)
            lp1 = diag(-1*ones(n(d)-1,1), [-1])+diag(2*ones(n(d), 1), [0])+diag(-1*ones(n(d)-1,1), [1]);
            Mp1 = diag((1/6)*ones(n(1)-1,1), [-1])+diag((4/6)*ones(n(1), 1), [0])+diag((1/6)*ones(n(1)-1,1), [1]);                                                    
            if (bound1(d)==1) lp1(1,1)=1; end;
            if (bound2(d)==1) lp1(n(d),n(d))=1; end;
            h = 1/(n(d)+1);
            if ((bound1(d)==1)||(bound2(d)==1)) h = 1/n(d); end;
            if ((bound1(d)==1)&&(bound2(d)==1)) h = 1/(n(d)-1); end;
            Mlp1 = diag(tt_a{d}(2:n(d),k), [-1])+diag(tt_a{d}(2:n(d),k), [1]);
            Mlp1 = Mlp1 + diag(tt_a{d}(1:n(d),k)+[tt_a{d}(2:n(d),k)' tt_a{d}(n(d),k)]', [0])*0.5;
        if (q==d)            
            cur_tt{d}(:,:,k) = Mlp1.*(lp1/h^2);
        else
            cur_tt{d}(:,:,k) = Mlp1.*Mp1;
%             cur_tt{d}(:,:,k) = diag(tt_a{d}(1:n(d),k)+[tt_a{d}(2:n(d),k)' tt_a{d}(n(d),k)]', [0])*0.5;
        end;
    end;    
    
    if (q==1)
        ttm = cur_tt;
    else
        ttm = ttm_add(ttm, cur_tt);
        ttm = tt_mat_compr(ttm, eps);
    end;
end;

end