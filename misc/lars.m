function [y]=lars(p,ind,mat,rhs,n,d,smin,h,i0,j0)
%The dimension is p^2 (thus w
k0=p*p;
x=zeros(k0,1);
pp=2.^(0:d-1);
for i=1:k0
   ix=ind(d*(i-1)+1:d*i);
   ix=dot((ix-1),pp); 
   x(i)=smin+h*ix;
end
mm=x(1)*mat{1};
for i=2:numel(mat)
  mm=mm+x(i)*mat{i};
end
   sol=mm \ rhs;
   sol=reshape(sol,n,n);
y=sol(i0,j0);
return
end