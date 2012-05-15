% test code for tt_qconv()
%
% April 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

% a,b,eps are parameters
d=10;
a=-1;
b=1;
eps=1E-10;
%

n=2^d;
h=(b-a)/n;
mesh=((1-n:n)'-1/2)*h;
x=exp(-1000*(mesh-0.1).^2);

mesh=((1-n/2:n/2)'-1/2)*h;
y=zeros(n,1);
y(n/2-n/8:n/2+n/8-1)=ones(n/4,1);
y=exp(-1000*(mesh-0.5).^2);


tt_x=full_to_tt(reshape(x,2*ones(1,d+1)),eps);
tt_y=full_to_tt(reshape(y,2*ones(1,d)),eps);

tt_z=tt_qconv(tt_x,tt_y,d);
z=h*tt_qtofull(tt_z,1);

f1=figure; hold on;
set(0,'CurrentFigure',f1)
plot(mesh,x(n/2:n*3/2-1),'r');
f2=figure; hold on;
set(0,'CurrentFigure',f2)
plot(mesh,y,'b');
f3=figure; hold on;
set(0,'CurrentFigure',f3)
plot(mesh,z,'g');