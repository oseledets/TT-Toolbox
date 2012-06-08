function [w]=wtt(v,wtt_transform)
%Forward WTT transform with previously computed filters
%   W=WTT(V,WTT_TRANSFORM) computes the forward WTT tranformation of a 
%   given tensor with filters, stored in the WTT_TRANSFORM structure.
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
u=wtt_transform.u;
r=wtt_transform.rks;
sz=wtt_transform.sz;

w=wtt_loc(v,u,r,sz);

return
end

function [w]=wtt_loc(v,u,r,sz)
if ( nargin == 1 )
  sz=size(v); %d=numel(sz);
else
  %d=numel(sz);
end
  N=numel(v);
  v=reshape(v,[r(1)*sz(1),N/(r(1)*sz(1))]);
if ( numel(u) == 1 )
    w=u{1}'*v;
    w=reshape(w,[r(1),N/(r(1))]);
else
  %Simplest one is recursion
  w=u{1}'*v;
  w0=w(1:r(2),:);
  unew=u(2:numel(u)); rnew=r(2:numel(r)); sznew=[sz(2),sz(3:numel(sz))];
  w0=wtt_loc(w0,unew,rnew,sznew);
  w(1:r(2),:)=w0;
  w=reshape(w,[r(1),N/r(1)]);
end
return
end
