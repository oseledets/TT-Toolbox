% test code for tt_qconv_p()
%
% April 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

% w,a,b,eps are parameters
d=12;
w=4;
a=0;
b=1;
r=10;
eps=1E-15;
%

n=2^d;
h=(b-a)/n;
mesh=((1:n)'-1/2)*h;

% x=sin(2*pi*w*mesh);
% %x=zeros(n,1);
% %x(n/2:n*3/4-1)=ones(n/4,1);
% 
% y=ones(n,1);
% y=(1:n)';
% y(n/4+1:n*3/4)=zeros(n/2,1);
% 
% tt_x=full_to_tt(reshape(x,2*ones(1,d)),eps);
% x=tt_qtofull(tt_x,1);
% tt_y=full_to_tt(reshape(y,2*ones(1,d)),eps);


tt_x=tt_rand(2,d,r);
tt_x=tt_compr2(tt_x,eps);
x=tt_qtofull(tt_x,1);
tt_y=tt_rand(2,d,r);
tt_y=tt_compr2(tt_y,eps);
y=tt_qtofull(tt_y,1);



tic;
tt_z=tt_qconv_p(tt_x,tt_y);
toc;
z=tt_qtofull(tt_z,1);
zzz=z;

tic;
tt=tt_qshiftstack_p(d);
tt=tt_qreshape(tt,3,[4*ones(d,1),2*ones(d,1)]);
%toc;

%tic;
tt=tt_mv(tt,tt_x);
%toc;


%tic;
tt=tt_qreshape(tt,1,2*ones(d,2));
tt_C=tt;
tt2_C=tt_matrix(tt_C);
%toc;

tt2_y=tt_tensor(tt_y);
%tt2_z=round(tt2_C,1e-12)*round(tt2_y,1e-12); tt2_z=round(tt2_z,1e-12);

%tic;
tt2_z=tt2_C*tt2_y;
toc;

z=full(tt2_z);
z=z(:);
ERRx=norm(zzz-z)/norm(z)




f1=figure; hold on;
set(0,'CurrentFigure',f1)
plot(mesh,x,'r');
f2=figure; hold on;
set(0,'CurrentFigure',f2)
plot(mesh,y,'b');
f3=figure; hold on;
set(0,'CurrentFigure',f3)
plot(mesh,z,'g');

