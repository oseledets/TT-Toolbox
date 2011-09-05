%Simple function to test the cross method
d=10;
%elem_fun=@(x) sum(x); %Just sum of everything 
p=0:d-1; p = 2.^p; 
a=-5;b=5;
n=2^d;
h=(b-a)/(n-1);
%mv=@(x) x.^3;
%elem_fun=@(x) 1.0./(dot((x-1),p)+1e-3); %Just sum of everything 
%elem_fun=@(x) mv(1e-12+dot(x-1,p)*h);

%Compare functions of TT-tensors
x=tt_x(d,2); x=tt_tensor(x); 
e=tt_ones(d,2);e=tt_tensor(e);
x=a*e+h*x; x=round(x,1e-13);
rs=x;
%fun=@(x) 1.0./sqrt(x);
%fun=@(x) 1.0./x;
fun=@(x) exp(-(x).^4) + 1;
elem_fun = @(ind) fun(rs(ind));
%elem_fun=@(ind) rs(ind);
%elem_fun=@(x) sqrt(x(1))+x(2);

eps=1e-6;
%y=tt_rc2(2*d,2,elem_fun,1e-12);
y=tt_rc(d,2,elem_fun,1e-6,'nswp',40,'change_dir_on',false);

z=funcrs2(rs,fun,1e-12,rs,20);


