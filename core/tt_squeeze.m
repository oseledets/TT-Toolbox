function [res]=tt_squeeze(tt)
%[RES]=TT_SQUEEZE(TT)
%Removes singleton dimensions in TT and returns reduced tensor RES
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
res=cell(d,1);
pos=1;
for k=1:d 
  if (size(tt{k},1) == 1 ) %Singleton dimension
     if ( k == 1 ) %Singleton dimension is the first one
       vec=tt{k}; r=size(vec,2); vec=reshape(vec,[1,r]);
       tt{k+1}=ten_conv(tt{k+1},2,vec'); 
       n=size(tt{k+1},1); r=size(tt{k+1},3);
       tt{k+1}=reshape(tt{k+1},[n,r]);
       
     elseif ( k == d ) %Singleton dimension is the last one
        vec=tt{k}; r=size(vec,2); vec=reshape(vec,[1,r]);
        res{pos}=ten_conv(res{pos},3,vec);
        n=size(res{pos},1); r=size(res{pos},2);
        res{pos}=reshape(res{pos},[n,r]);
     else
         mat=tt{k}; r1=size(mat,2); r2=size(mat,3);
        mat=reshape(mat,[r1,r2]);
        tt{k+1}=ten_conv(tt{k+1},2,mat');
        r3=size(tt{k+1},3); n=size(tt{k+1},1);
        tt{k+1}=reshape(tt{k+1},[n,r1,r3]);
     end
  else
      res{pos}=tt{k};
      pos=pos+1;
  end
end
pos=pos-1;
res=res(1:pos);
return
end
