%Simple function to test the cross method
d=4;
%elem_fun=@(x) sum(x); %Just sum of everything 
p=0:d-1; p = 2.^p; 
a=0;b=1;
n=2^d;
h=(b-a)/(n-1);
mv=@(x) x.^2;
%elem_fun=@(x) 1.0./(dot((x-1),p)+1e-3); %Just sum of everything 
elem_fun=@(x) mv(dot(x-1,p)*h);
eps=1e-6;
y=tt_rc(d,2,elem_fun,eps);



