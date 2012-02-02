function [D]=lagr_diff(x)
% function [D]=lagr_diff(x)
% 
% Generates the differentiation matrix in the basis of the Lagrange
% interpolating polynomials on grid X, so that D(i,j) = L'_j(x_i)

n = max(size(x));
D = zeros(n,n);

denom = ones(n,1);
for i=1:n % grid points in result
    for k=1:n % inner product
        if (k~=i)
            denom(i)=denom(i)*(x(i)-x(k));
        end;
    end;
end;

for j=1:n % grid points in result
    for i=1:n % polynomial index
        curD=0;
        for p=1:n % sum
            if (p~=i)
                curnom = 1;
                for k=1:n % product
                    if (k~=i)&&(k~=p)
                        curnom=curnom*(x(j)-x(k));
                    end;
                end;
                curD=curD+curnom;
            end;
        end;
        D(j, i)=curD/denom(i);
    end;
end;

end