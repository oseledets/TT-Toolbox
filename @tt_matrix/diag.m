function [tt]=diag(tm)
%[TM]=DIAG(TT)
%Extracts diagonal of the TT-matrix TM
tt=tt_tensor;
n=tm.n;
tt1=tm.tt;
r=tt1.r;
core=tt1.core;
ps=tt1.ps;
d=tt1.d;
ps1=cumsum([1;n.*r(1:d).*r(2:d+1)]);
core1=zeros(ps1(d+1)-1,1);
for i=1:d
   cur_core=core(ps(i):ps(i+1)-1);cur_core=reshape(cur_core,[r(i),n(i),n(i),r(i+1)]);
   cur_core1=zeros(r(i),n(i),r(i+1));
   for j=1:n(i)
      cur_core1(:,j,:)=cur_core(:,j,j,:);
   end
   core1(ps1(i):ps1(i+1)-1)=cur_core1(:);
end
tt.r=r;
tt.n=n;
tt.core=core1;
tt.ps=ps1;
tt.d=d;
return
end