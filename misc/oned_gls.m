function [A,Aadd,b,badd] = oned_gls(k, v, s, f, N, a, b)

h = (b-a)/(N+1);
xint = (a+h:h:b-h)';
% xhalf = (a+h/2:h:b-h/2)';
xintb = (a:h:b)';

% kvec = k(xhalf);
kvec = (k(xintb(1:end-1))+k(xintb(2:end)))*0.5;
kvec(isinf(kvec))=0;
kvec(isnan(kvec))=0;

% vvec = v(xhalf);
vvec = (v(xintb(1:end-1))+v(xintb(2:end)))*0.5;
vvec(isinf(vvec))=0;
vvec(isnan(vvec))=0;

% svec = s(xhalf);
svec = (s(xintb(1:end-1))+s(xintb(2:end)))*0.5;
svec(isinf(svec))=0;
svec(isnan(svec))=0;

% fvec = f(xhalf);
fvec = (f(xintb(1:end-1))+f(xintb(2:end)))*0.5;
fvec(isinf(fvec))=0;
fvec(isnan(fvec))=0;

M = spdiags(ones(N,1)*[1/6, 4/6, 1/6]*h, [-1,0,1], N, N);
S1 = spdiags(ones(N,1)*[-0.5, 1, 0.5], [-1,0,1], N, N);
S2 = spdiags(ones(N,1)*[-1, 2, -1]/h, [-1,0,1], N, N);

% tau-vector
%  tau = 4*kvec/h^2 + 2*abs(vvec)/h + abs(svec);
%  tau = 1./tau;
tau = ones(N+1,1);

% Prepare averaged diagonals
% Mass-style - need for k, s, v.*k, v.*v
kdiag = zeros(N,3);
kdiag(1:N-1,1) = kvec(2:N);
kdiag(1:N,2) = (kvec(1:N)+kvec(2:N+1))*0.5;
kdiag(2:N,3) = kvec(2:N);
kdiag = spdiags(kdiag, [-1, 0, 1], N, N);

sdiag = zeros(N,3);
sdiag(1:N-1,1) = svec(2:N);
sdiag(1:N,2) = (svec(1:N)+svec(2:N+1))*0.5;
sdiag(2:N,3) = svec(2:N);
sdiag = spdiags(sdiag, [-1, 0, 1], N, N);

% vkdiag = zeros(N,3);
% vkdiag(1:N-1,1) = vvec(2:N).*kvec(2:N).*tau(2:N);
% vkdiag(1:N,2) = (vvec(1:N).*kvec(1:N).*tau(1:N)+vvec(2:N+1).*kvec(2:N+1).*tau(2:N+1))*0.5;
% vkdiag(2:N,3) = vvec(2:N).*kvec(2:N).*tau(2:N);
% vkdiag = spdiags(vkdiag, [-1, 0, 1], N, N);

vvdiag = zeros(N,3);
vvdiag(1:N-1,1) = vvec(2:N).*vvec(2:N).*tau(2:N);
vvdiag(1:N,2) = (vvec(1:N).*vvec(1:N).*tau(1:N)+vvec(2:N+1).*vvec(2:N+1).*tau(2:N+1))*0.5;
vvdiag(2:N,3) = vvec(2:N).*vvec(2:N).*tau(2:N);
vvdiag = spdiags(vvdiag, [-1, 0, 1], N, N);


% Grad-style - need for v, v.*s
vdiag = zeros(N,3);
vdiag(1:N-1,1) = vvec(2:N);
vdiag(1:N,2) = (vvec(1:N)-vvec(2:N+1))*0.5;
vdiag(2:N,3) = vvec(2:N);
vdiag = spdiags(vdiag, [-1, 0, 1], N, N);

vsdiag = zeros(N,3);
vsdiag(1:N-1,1) = vvec(2:N).*svec(2:N).*tau(2:N);
vsdiag(1:N,2) = (vvec(1:N).*svec(1:N).*tau(1:N)-vvec(2:N+1).*svec(2:N+1).*tau(2:N+1))*0.5;
vsdiag(2:N,3) = vvec(2:N).*svec(2:N).*tau(2:N);
vsdiag = spdiags(vsdiag, [-1, 0, 1], N, N);

Pe = abs(v(xint))*h./(2*k(xint));
tau = h./(2*abs(v(xint))).*(coth(Pe)-1./Pe);
%  Pe = abs(vvec)*h./(2*kvec);
%  tau = h./(2*abs(vvec)).*(coth(Pe)-1./Pe);
tau(isnan(tau))=0;
tau(isinf(tau))=0;
%  tau = 4*k(xint)/h^2 + 2*abs(v(xint))/h + abs(s(xint));
%  tau = 4*kvec/h^2 + 2*abs(vvec)/h + abs(svec);
%  tau = 1./tau;
%  tau = (tau(1:N)+tau(2:N+1))*0.5;
%  tau = 1./tau
tau = spdiags(tau, 0, N, N);

% Simple rhs
b = (fvec(1:N)+fvec(2:N+1))*0.5*h;
% Add-rhs
%  badd = (vvec(1:N).*fvec(1:N).*tau(1:N)-vvec(2:N+1).*fvec(2:N+1).*tau(2:N+1))*0.5;
badd = tau*(vvec(1:N).*fvec(1:N)-vvec(2:N+1).*fvec(2:N+1))*0.5;

A = kdiag.*S2 + vdiag.*S1 + sdiag.*M;
% Stab. matrix
Aadd = tau*((vvdiag).*S2 + vsdiag.*(S1'));

end