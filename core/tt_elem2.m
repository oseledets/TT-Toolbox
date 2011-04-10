function [res]=tt_elem2(tt,ind)
%[RES]=TT_ELEM2(TT,IND)
%Returns the element of TT tensor in position given by multiindex IND
%If IND is short of length l < d it is assumed it goes from d-l+1 to d
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
mat=tt{d};
len=max(size(ind));
v0=mat(ind(len),:); v0=v0';
ps=len-1;
if ( len < d )
 for k=d-1:-1:d-len+1
  core=tt{k};
  r2=size(core,2); r3=size(core,3);
  core=core(ind(ps),:,:);
  core=reshape(core,[r2,r3]);
  %core=squeeze(core); %To avoid rank-1 terms
  %v0=squeeze(core(ind(ps),:,:))*v0;
  v0=core*v0;
  %v0=v0';
  ps=ps-1;
 end
 res=tt{1}*v0;
else
  res=tt_elem(tt,ind);    
end
return
end
