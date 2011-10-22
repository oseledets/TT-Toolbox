function [D,x]=spectral_chebdiff(N, a, b)
%[D,x]=spectral_chebdiff(N,a,b)
%Computes the Spectral Chebyshev discretization
%of on the interval [a,b], together with the 
%interpolation nodes
x = a+(b-a)*0.5*(cos(pi*(0:N)/N)+1)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = (X-X')/(b-a);
D = (c*(1./c)')./(dX+eye(N+1));
D = D-diag(sum(D, 2));
D = D/(b-a);

end