function [tt]=tt_squeeze(tt)
%Removes singleton dimensions from the TT-tensor
%   [TT]=TT_SQUEEZE(TT) Removes singleton dimensions from the TT-tensor.
%   Uses TT1.0 format. Please avoid its usage: it will be removed in
%   future releases. Use reshape() and squeeze() from the object-oriented 
%   version
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

% 
% If a singleton dimension is in the middle of a train - everything breakes
% Use tt_reshape instead!

d=size(tt,1);
k=1;
while ( k <= d )
  if (size(tt{k},1) == 1 ) %Singleton dimension
     if ( k == 1 ) %Singleton dimension is the first one
       vec=tt{k}; r=size(vec,2); vec=reshape(vec,[1,r]);
       tt{k+1}=ten_conv(tt{k+1},2,vec'); 
       n=size(tt{k+1},1); r=size(tt{k+1},3);
       tt{k+1}=reshape(tt{k+1},[n,r]);
       tt{k}=[];
     else
         %Delete the dimension k (and incorporate into the
         %k-1)        
         mat=tt{k}; r1=size(mat,2); r2=size(mat,3);
        mat=reshape(mat,[r1,r2]);
        tt{k-1}=ten_conv(tt{k-1},3,mat);
        tt(k)=[];
        d=d-1;
     end
  else
     k=k+1;    
  end
end
return
end
