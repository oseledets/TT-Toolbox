%Simple function to test the cross method
d=5;
n=2;
elem_fun=@(x) sum(x); %Just sum of everything 
eps=1e-6;
y=tt_rc(d,n,elem_fun,eps);



