function [M]=poly_bilin(c, wfun, D1, D2)

n = size(c,1)-1;
[nodes, cw]=fclencurt(n+1, -1, 1);
w = wfun(nodes);

M = zeros(n+1,n+1);
for j=0:n
    for i=0:n
        p1 = build_poly(i, nodes, c);
        p2 = build_poly(j, nodes, c);
        p1 = D1*p1;
        p2 = D2*p2;
        M(i+1,j+1) = sum(p1.*p2.*w.*cw);
    end;
end;

end