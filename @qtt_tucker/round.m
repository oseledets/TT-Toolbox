function [tt]=round(tt,varargin)
%[QTT_TUCKER]=ROUND(QTT_TUCKER,EPS)
%The rounding procedure for the QTT-Tucker class
d=tt.dphys;
core=tt.core;
tuck=tt.tuck;
eps=varargin{1};
for i=1:d
   [tuck{i},rm]=qr(tuck{i},'lr');
   core{i}=ten_conv(core{i},2,rm.');
end
core=round(core,eps); %Round the core --- we know the result comes
%with lr orthogonality?
rtt=rank(core); 
n=size(core);
for i=1:d
   cr=reshape(core{i},[rtt(i),n(i),rtt(i+1)]);
   cr=permute(cr,[2,1,3]); cr=reshape(cr,n(i),rtt(i)*rtt(i+1));
   [u,s,v]=svd(cr,'econ');
   s=diag(s);
   r=my_chop2(s,norm(s)*eps);
   u=u(:,1:r); s=s(1:r); v=v(:,1:r);
   tuck{i}=tuck{i}*(u*diag(s)); 
   tuck{i}=round(tuck{i},eps);
   [tuck{i},rm]=qr(tuck{i},'lr');
   cr=rm*v';
   cr=reshape(cr,[r,rtt(i),rtt(i+1)]);
   core{i}=permute(cr,[2,1,3]);
   %The result is lr-orthogonal, thus we need only to 
end
   
return
end
