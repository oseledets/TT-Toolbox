d=10; n=2^d;
a=0; b=1; h=(b-a)/(n-1);
x=a+(0:n-1)*h; true_x=tt_tensor(reshape(x,2*ones(1,d)),1e-8);
x_p=x;
x_p(1:n/2)=randn(n/2,1);
xx=tt_tensor(reshape(x_p,2*ones(1,d)),1e-12); 
%Generate weight
w=ones(1,n);
w(1:n/2)=0;
ww=tt_tensor(reshape(w,2*ones(1,d)),1e-12);
%wround(x,ww,
norm(diag(ww)*(xx-true_x))
xxx=wround(xx,ww,1e-8,[],5);
plot(full(xxx,[n,1]))