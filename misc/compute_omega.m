function [w1,w2]=compute_omega(a,c,N)
%[w1,w2]=compute_omega(a,c,N);
%Finds first N roots roots to the two equations
%c==w*tan(a*w) (w1)
%w + c*tan(a*w) == 0 (w2) 
%in an interval [0,q]

%Compute bracketing intervals: 
%pi/(2*a)*(2*N-1) = q
%2*N-1 = (2*q*a)/(pi)
%N=floor(q*a/pi+0.5);
w=pi/(a*2)*(1:2:2*N+3);
ww1=[0,w(1:N)];
ww2=w(1:N+1);
dt=1e-13;
tol=1e-12;
w1=[];
w2=[];
for i=1:N
  w1(i)=illinois(@(x) x*tan(a*x) - c, ww1(i)+dt, ww1(i+1)-dt,tol);
  w2(i)=illinois(@(x) x + c * tan(a*x), ww2(i)+dt, ww2(i+1)-dt,tol);
end

return
end