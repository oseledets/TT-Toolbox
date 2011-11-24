function [tt]=ttm(tt,k,mat)
%[tt]=ttm(tt,k,mat)
%Multiplies TT-tensor over the K-th mode with the matrix MAT

cr=tt.core; 
if (size(cr,1) == 1)
  cr=cr.'; %This procedure needs "column" cr
end
ps=tt.ps;
d=tt.d;
n=tt.n;
r=tt.r;
m=size(mat,2);
crx=cr(ps(k):ps(k+1)-1);
crx=reshape(crx,[r(k),n(k),r(k+1)]);
crx=permute(crx,[1,3,2]); 
crx=reshape(crx,[r(k)*r(k+1),n(k)]);
crx=crx*mat;
%crs is r(k)xr(k+1)xm
crx=reshape(crx,[r(k),r(k+1),m]);
crx=permute(crx,[1,3,2]);
cr=[cr(1:ps(k)-1);crx(:);cr(ps(k+1):end)];
n(k)=m;
tt.n=n;
tt.r=r;
tt.ps=cumsum([1;n.*r(1:d).*r(2:d+1)]);
tt.core=cr;
return
end