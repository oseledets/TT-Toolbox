function [tt]=canm_to_qtt(can,rmax,eps)
%[TT]=CANM_TO_QTT(CAN,RMAX,EPS)
% Converts canonical representation CAN (Cell array of factors)
% of a matrix
% to QTT representation with accuracy parameter EPS
% Mode sizes should be powers of 2
%
%
% TT Toolbox 1.0
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
tt=tt_zeros(sum(d),4);
tt_arr=cell(d1,1);
for i=1:R
  for s=1:d1
    %tt_arr{s}=vec_to_nd(can{s}(:,i),eps,d(s));
    tt_arr{s}=tt_mat_to_vec(full_to_nd(can{s}(:,:,i),eps,d(s))); 
 end
  mat=tt_mkron(tt_arr);
  tt=tt_add(tt,mat);
  if ( tt_erank(tt) > rmax )
     tt=tt_compr2(tt,eps);
  end
end
tt=tt_compr2(tt,eps);
return
