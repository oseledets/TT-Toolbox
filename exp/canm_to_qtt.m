function [tt]=canm_to_qtt(can,rmax,eps)
%Converts a canonical representation to QTT-format
%   [TT]=CANM_TO_QTT(CAN,RMAX,EPS)
%   Converts canonical representation CAN (Cell array of factors)
%   of a matrix
%   to QTT representation with accuracy parameter EPS
%   Mode sizes should be powers of 2
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

if (size(can,1) == 1 )
  can=can';
end
d1=size(can,1);
R=size(can{1},3);
d=zeros(d1,1);
for i=1:d1
  w1=size(can{i},1);
    w=log2(w1); w=round(w);
  if (2^w ~= w1 )
     fprintf('can_to_qtt error: mode %d size = %d not a power of two \n',i,w1);
     return
  end
  d(i)=w;
end
%For each factor generate its qtt-representation
tt=[];
tt_arr=cell(d1,1);
for i=1:R
  for s=1:d1
    %tt_arr{s}=vec_to_nd(can{s}(:,i),eps,d(s));
    can_cur=can{s}(:,:,i); 
    can_cur=reshape(can_cur,2*ones(1,2*d(s)));
    tt_arr{s}=tt_tensor(tt_matrix(can_cur,eps));
  end
  mat=tt_mkron(tt_arr);
  tt=tt+mat;
  if ( erank(tt) > rmax )
     tt=round(tt,eps);
  end
end
tt=round(tt,eps);
tt=tt_matrix(tt);
return
