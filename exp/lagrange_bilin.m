function [M]=lagrange_bilin(xb, cw, wfun, D1, D2)

N = max(size(xb));
if (isa(wfun, 'double'))&&(max(size(wfun))==N)
    W=wfun;
else
    W=wfun(xb);
end;

M = zeros(N,N);
for j=1:N
    for i=1:N
        M(i,j)=sum(W.*cw.*D1(:,i).*D2(:,j));
    end;
end;

end