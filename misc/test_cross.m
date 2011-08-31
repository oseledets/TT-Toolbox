%Simple function to test the cross method
d=10;
%elem_fun=@(x) sum(x); %Just sum of everything 
p=0:d-1; p = 2.^p; 
a=-5;b=5;
n=2^d;
h=(b-a)/(n-1);
%mv=@(x) x.^3;
mv=@(x) sqrt(x)+abs(x);
%elem_fun=@(x) 1.0./(dot((x-1),p)+1e-3); %Just sum of everything 
%elem_fun=@(x) mv(1e-12+dot(x-1,p)*h);

%Compare functions of TT-tensors
x=tt_x(d,2); x=tt_tensor(x); 
e=tt_ones(d,2);e=tt_tensor(e);
x=a*e+h*x; x=round(x,1e-13);
x1=kron(e,x); x2=kron(x,e);
rs=x1.^2+x2.^2; rs=round(rs,1e-13);
%fun=@(x) 1.0./sqrt(x);
%fun=@(x) 1.0./x;
fun=@(x) exp(-x.^4);
%fun = @(x) x.^2;
rs=x;
%elem_fun = @(ind) fun(rs(ind));
%elem_fun=@(ind) rs(ind);
%elem_fun=@(x) sqrt(x(1))+x(2);
elem_fun = @(ind) fun(rs(ind));
eps=1e-6;
D=ndims(rs);
%y=tt_rc2(2*d,2,elem_fun,1e-12);
y=tt_rc(D,2,elem_fun,1e-6,'nswp',40);
%v=tt_rand(size(y),ndims(y),2);
%y=y+v; y=round(y,1e-12);
%y1=tt_rc(D,2,elem_fun,1e-8,'nswp',40,'x0',y);

z=funcrs2(rs,fun,1e-6,rs,20);
z=round(z,1e-12);
z1=tt_rc(D,2,elem_fun,1e-8,'nswp',40,'x0',z);


