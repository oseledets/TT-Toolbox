function [res]=tt_mkron(tt_arr)
%[RES]=TT_MKRON(TT_ARR)
% Computes Kronecker product of many TT tensors
% TT_ARR is a cell array of TT tensors,
% and the result RES is the Kronecker product of them
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
if ( size(tt_arr,1) == 1 ) 
  tt_arr=tt_arr';
end
N=size(tt_arr,1);
dsz=zeros(N,1);
for k=1:N
    if ( iscell(tt_arr{k}) ) 
      dsz(k) = size(tt_arr{k},1);
    else
      dsz(k) = 1;     
    end
end
d=sum(dsz(1:N));
res=cell(d,1);
pos=1;
for k=1:N
  for i=1:dsz(k)
      if ( ~iscell(tt_arr{k}))
        res{pos}=tt_arr{k};
      else
     res{pos}=tt_arr{k}{i};
      end
     pos=pos+1;
  end
end
pos=dsz(1)+1;
for k=2:N
    fx=res{pos};
    ncur=size(fx,1); r=size(fx,2);
    res{pos}=reshape(fx,[ncur,1,r]);
    pos = pos + dsz(k);
end

return
end
