function [newtt, rv] = tt_blqr(tt)
%[NEWTT,RV]=TT_BLQR(TT,EPS)
%Reorthogonalizes the colunms in the block
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
d = size(tt,1);
newtt=cell(d,1);

%First, estimate the logarithm of the squared norm, and then --- the norm
%of each core :-)
nrm_full=tt_dot2(tt,tt); 
if ( nrm_full < log(1e-200) ) %Something really small
  n=size(tt{1},1);
  newtt{1}=zeros(n,0);
  n=size(tt{d},1);
  newtt{d}=zeros(n,0);
  for i=2:d-1
      n=size(tt{i},1);
      newtt{i}=zeros(n,0,0);
  end
  rv = [];
  return
end
%nrmf=zeros(d,1); %Place to store the norms of QR-factors
%first, we orthogonalize the tensor from left to right until the last mode

[q,rv]=qr(tt{1},0); 

tt{1}=q;
for i=2:d-1
    tt{i} = ten_conv(tt{i},2,rv');
    ncur=size(tt{i},1);
    r2=size(tt{i},2);
    r3=size(tt{i},3);
    core=reshape(tt{i},[ncur*r2,r3]);
    [tt{i},rv]=qr(core,0); 
% %    rv=rv./nrm1;
    rnew=min(r3,ncur*r2);
% % tt{i}=reshape(tt{i}.*nrm1,[ncur,r2,rnew]);
    tt{i}=reshape(tt{i},[ncur,r2,rnew]);
end

tt{d} = ten_conv(tt{d},2,rv');
 [core,rv]=qr(tt{d}',0); 
 tt{d}=core';
% [core,rv]=qr(tt{d},0); 
% tt{d}=core;
% rv=rv';

newtt=tt;

return
end
