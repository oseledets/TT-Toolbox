function [w]=iwtt(v,wtt_transform)
u=wtt_transform.u;
r=wtt_transform.rks;
sz=wtt_transform.sz;
w=iwtt_loc(v,u,r,sz);
return
end
function [w]=iwtt_loc(v,u,r,sz)
%[W]=IWTT(V,U,R,SZ)
%Applies Inverse WTT transform to vector V with linear filters U 
%and ranks R
%SZ is optional
if ( nargin == 1 )
  sz=size(v);% d=numel(sz);
else
 % d=numel(sz);
end
%Apply one transform back
  N=numel(v);
  v=reshape(v,[r(1)*sz(1),N/(r(1)*sz(1))]);
if ( numel(u) == 1 )
   w=u{1}*v;
    w=reshape(w,[r(1),N/(r(1))]);
else
  %Simplest one is recursion
  w=v;
  w0=v(1:r(2),:);
  unew=u(2:numel(u)); rnew=r(2:numel(r)); sznew=[sz(2),sz(3:numel(sz))];
  w0=iwtt_loc(w0,unew,rnew,sznew);
  w(1:r(2),:)=w0;
  m=size(u{1},1);
  w=reshape(w,[m,N/m]);
  w=u{1}*w;
  w=reshape(w,[r(1),N/r(1)]);
end
  return
end
