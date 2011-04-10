function [a] = full(tt)
%[A]=FULL(TT)

d=tt.d; n=tt.n; ps=tt.ps; core=tt.core; r=tt.r;
a=1;
for i=1:d
  cr=core(ps(i):ps(i+1)-1);
  cr=reshape(cr,[r(i),n(i)*r(i+1)]);
  a=reshape(a,[numel(a)/r(i),r(i)]);
  a=a*cr;
end
a=reshape(a,n');

return;
