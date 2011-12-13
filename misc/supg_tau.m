function [tau] = supg_tau(v, k, h)
%  function [tau] = supg_tau(v, k, h)
%  Computes the SUPG coefficients tau = h/(2*|v|)*(coth(Pe)-1./Pe)
%  accurately, using the Laurent expansion if |v|*h<tol

tol = 1e-4;

v = v(:);
k = k(:);

n = max(numel(v), numel(k));
if (numel(k)==1)
    k = k*ones(n,1);
end;
if (numel(v)==1)
    v = v*ones(n,1);
end;

Pe = abs(v)*h./(2*k);
ctz = find(abs(v*h)<tol); % small elements are to approximate via Laurent series
not_ctz = find(abs(v*h)>=tol); % direct comput. is possible

tau = zeros(n,1);
tau(not_ctz) = h*0.5*(coth(Pe(not_ctz))-1./Pe(not_ctz))./abs(v(not_ctz));
tau(ctz) = h^2./(12*k(ctz)) - h^4*(v(ctz).^2)./(720*k(ctz));

end