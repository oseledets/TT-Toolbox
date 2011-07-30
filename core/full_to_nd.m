function [res]=full_to_nd(mat,eps,d1,d2) 
%[RES]=FULL_TO_ND(MAT,EPS,D)
%Converts 2^D1 x 2^D2 matrix into QTT format with accuracy EPS 
% TT Toolbox 1.0, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
if (nargin<3)||(isempty(d1))
    d1 = log2(size(mat,1));
end;
if (nargin<4)||(isempty(d2))
    d2 = log2(size(mat,2));
    if (nargin>2)&&(~isempty(d1))
        d2 = d1;
    end;
end;
d = max(d1,d2);
sz=[2*ones(1,d1), ones(1,max(d2-d1,0)), 2*ones(1,d2), ones(1,max(d1-d2,0))];
res=reshape(mat,sz); %(i_1,...,i_d0,j_1,...,j_d0)
prm=zeros(1,d*2);
for i=1:d
  prm(2*i-1)=i;
  prm(2*i)=d+i;
end
res=permute(res,prm);
sz1 = prod(reshape(sz(prm), 2, d), 1);
% sz1=4*ones(1,d);
res=reshape(res,sz1);
res=full_to_tt(res,eps);
res=tt_vec_to_mat(res,sz(1:d)',sz(d+1:2*d)');
return
end
