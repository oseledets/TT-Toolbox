function [tt] = multigauss(P, n, a, eps)
% function [tt] = multigauss(P, n, a, eps)
% Generates a generalized multidimensional Gaussian function
% exp(-x'Px) on a uniform grid in the TT format
% n is the grid sizes,
% a is the domain size (domain is [-a; a]^d),
% eps - approximation tolerance for TT truncations

d = size(P,1);

if (numel(n)==1)
    n = n*ones(1,d);
end;
if (numel(a)==1)
    a = a*ones(1,d);
end;

tt = tt_tensor(tt_ones(d, n));

for i=1:d
    for j=1:d        
        if (i~=j)
            h1 = 2*a(i)/n(i); h2 = 2*a(j)/n(j);
            x1 = (-a(i):h1:a(i)-h1)'*ones(1,n(j));
            x2 = ones(n(i),1)*(-a(j):h2:a(j)-h2);
            
            exp2 = exp(-x1.*P(i,j).*x2);
            exp2 = full_to_tt(exp2, eps);
            expcr1 = exp2{1};
            r = size(expcr1, 2);
            expcr2 = exp2{2};
            ex1 = tt_tensor; ex1.n = n(i); ex1.d = 1; ex1.ps = cumsum([1; n(i)*r]);
            ex2 = tt_tensor; ex2.n = n(j); ex2.d = 1; ex2.ps = cumsum([1; n(j)*r]);
            if (i<j)
                expcr1 = reshape(expcr1, 1, size(expcr1,1), size(expcr1, 2));
                expcr2 = expcr2.';
                ex1.r = [1;r]; ex2.r = [r;1];
            else
                expcr1 = expcr1.';
                expcr2 = reshape(expcr2, 1, size(expcr2,1), size(expcr2, 2));
                ex1.r = [r;1]; ex2.r = [1;r];
            end;
            ex1.core = expcr1(:);
            ex2.core = expcr2(:);
        else
            h1 = 2*a(i)/n(i);
            x1 = (-a(i):h1:a(i)-h1)';
            exp2 = exp(-x1.*P(i,i).*x1);
        end;
        
        curtt = [];
        for k=1:d
            e1 = tt_tensor(tt_ones(1, n(k)));
            if (k<min(i,j))||(k>max(i,j))
                curtt = kron(curtt, e1);
            elseif (k==i)&&(k==j)
                curtt = kron(curtt, tt_tensor(exp2));
            elseif (k==i)&&(k~=j)
                curtt = kron(curtt, ex1);
            elseif (k==j)&&(k~=i)
                curtt = kron(curtt, ex2);
            elseif (k>min(i,j))&&(k<max(i,j))&&(i~=j)
                e1r = tt_tensor; e1r.d =1; e1r.n = n(k); e1r.r = [r;r];
                e1r.ps = cumsum([1; r^2*n(k)]); 
                e1cr = kron(eye(r), ones(n(k),1));
                e1cr = reshape(e1cr, n(k),r, r);
                e1cr = permute(e1cr, [2, 1, 3]);
                e1r.core = e1cr(:);
                curtt = kron(curtt, e1r);
            end;
        end;
        tt = tt.*curtt;
        tt = round(tt, eps);
        fprintf('=multigauss= term (%d,%d) is done. Erank: %g\n', i, j, erank(tt));
    end;
end;

end