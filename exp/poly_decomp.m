function [alpha]=poly_decomp(ffun, c, wfun)

n = size(c,1)-1;
alpha = zeros(n+1,1);
[nodes, cw]=fclencurt(n+1, -1, 1);
f = ffun(nodes); w = wfun(nodes);

for i=1:n+1
    p = build_poly(i-1, nodes, c);
    alpha(i)=sum(p.*f.*w.*cw);
end;

end