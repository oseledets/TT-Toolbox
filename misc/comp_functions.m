%Compute Karhunen-Loeve expansion for the function exp(-|x1-x2|  - |y1-y2|)

p=4; %The number of 1d parameters
a0=1; %Interval is [-a0,a0]
c0=1; %Inverse correlation length
[w1,w2]=compute_omega(a0,c0,p);
%Compute eigev
lam1=2*c0./(w1.^2+c0.^2);

lam2=2*c0./(w2.^2+c0.^2);


L=6; n=2^L-1; m=(n+1)/p; %Step size
a=-a0;
b=a0;
%The coefficient has to be specified in the middle points

t=(0.5:n+0.5)/(n+1);
x=a+(b-a)*t;

%Compute one-dimensional eigenfunctions on this grid
f=zeros(n+1,2*p); 
lam=zeros(1,2*p);
for i=1:p
f(:,2*i-1)=cos(w1(i)*x)./sqrt(a0+sin(2*w1(i)*a0)./(2*w1(i)));
f(:,2*i)  =sin(w2(i)*x)./sqrt(a0-sin(2*w2(i)*a0)./(2*w2(i)));
lam(2*i-1)=lam1(i);
lam(2*i)=lam2(i);
end

