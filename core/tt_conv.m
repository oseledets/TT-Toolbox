function [res] = tt_conv(tt,vec)
%[RES]=TT_CONV(TT,VEC)
% Computes A x VEC{1} x ... x VEC{D}
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
d=size(tt,1);
mat=cell(d,1);

v0=vec{1}'*tt{1};
vd=vec{d}'*tt{d};
for i=2:d-1
  core=tt{i};
  ncur=size(core,1);
  r2=size(core,2);
  r3=size(core,3);
  core=reshape(core,[ncur,r2*r3]);
  mat{i}=reshape(vec{i}'*core,[r2,r3]);
end

for i=2:d-1
  v0=v0*mat{i};  
end
res=v0*vd';
return
end
