%Simple function to test the cross method
d=20;
%elem_fun=@(x) sum(x); %Just sum of everything 
p=0:d-1; p = 2.^p; 
a=0;b=1;
n=2^d;
h=(b-a)/(n-1);
%mv=@(x) x.^3;
mv=@(x) sqrt(x)+abs(x);
%elem_fun=@(x) 1.0./(dot((x-1),p)+1e-3); %Just sum of everything 
%elem_fun=@(x) mv(1e-12+dot(x-1,p)*h);

%Compare functions of TT-tensors
x=tt_x(d,2); x=tt_tensor(x); 
e=tt_ones(d,2);e=tt_ones(e);
%elem_fun=@(x) sqrt(x(1))+x(2);

eps=1e-6;
y=tt_rc(d,2,elem_fun,eps);



