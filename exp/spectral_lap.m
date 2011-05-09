%Harmonic oscillator using Chebyshev discretization by Trefethen
K=4;
format long, format compact;
L = 20;% domain is [-L L], periodic
  cv=[]; q=1;
function [lp,x]=spectral_lap(N,L)
%[LP,X]=SPECTRAL_LAP(N,L)
%Creates a Chebyshev-type discretization of the Laplace operator on a Chebyshev grid
%in [-L,L], with periodic BC
%Also, the grid itself is returned
%Typicall, N can be small
h = 2*pi/N; x = h*(1:N); x = L*(x-pi)/pi;
    column = [-pi^2/(3*h^2)-1/6 ...
       -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
    lp = (pi/L)^2*toeplitz(column);  % 1d laplace
    lp=-lp;
return
