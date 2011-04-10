function [tt] = core(tt1,varargin)
%[TT]=CORE(TT1)
%[TT]=CORE(TT1,NUM)
%Converts TT-tensor TT1 to old-cell array format,
%optionally extract the required core
if (nargin==1)
  tt2=tt1.tt;
  n=tt1.n; m=tt1.m;
  tt=core(tt2);
  tt=tt_vec_to_mat(tt,n,m);
else
  i=varargin{1};
  tt2=tt1.tt; 
  tt=core(tt2,i);
  r=tt2.r;
  n=tt1.n; m=tt1.m;
  tt=reshape(tt,[r(i),n(i),m(i),r(i+1)]);
end
return
end