function [p]=build_poly(n, nodes, c)

p = 0*nodes;
for i=n:-1:0
%     p = p+c(n+1,i+1)*(nodes.^i);
    p = p+c(n+1,i+1)*cos(i*acos(nodes));
end;

end