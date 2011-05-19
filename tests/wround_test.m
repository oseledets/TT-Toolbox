d=10; n=2^d;
a=0; b=1; h=(b-a)/(n-1);
x=a+(0:n-1)*h; true_x=tt_tensor(reshape(x,2*ones(1,d)),1e-8);
e=tt_ones(d,2); e=tt_tensor(e);
x_p=x;
x_p(1:n/2)=randn(n/2,1);
xx=tt_tensor(reshape(x_p,2*ones(1,d)),1e-12); 
%Generate weight
w=ones(1,n);
w(1:n/2)=0;
ww=tt_tensor(reshape(w,2*ones(1,d)),1e-12);
%wround(x,ww,
norm(diag(ww)*(xx-true_x))
%xxx=true_x+1e-1*xx;
xxx=xx;
close all;
for i=1:10000

xxx1=round(ww.*xx+(e-ww).*xxx,1e-6,2);
if ( mod(i,400) == 0 )
fprintf('i=%d er=%3.2e \n',i,norm(xxx1-xxx)/norm(xxx));
plot(full(xxx,[n,1]));
drawnow;
end
xxx=xxx1;

end
%xxx=wround(xx,ww,1e-8,[],5);
plot(full(xxx,[n,1]))