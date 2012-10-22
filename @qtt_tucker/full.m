function [a] = full(tt, sizes)
%Converts a QTT-Tucker tensor the full tensor
%   [A]=FULL(QTT_TUCKER,SIZES) Converts QTT_Tucker tensor into a tensor
%   with mode sizes SIZES. The dimensions should be consistent.
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

n=size(tt.core);
r = tt.core.r;
a=core2cell(tt.core); 
fac=cell(tt.dphys,1);
for i=1:tt.dphys
   ti=tt.tuck{i}; ri=rank(ti,ti.d+1);
    ti=full(ti);
   fac{i}=reshape(ti,[numel(ti)/ri,ri]);
end
%And now convert tucker representation to the full one
for i=1:tt.dphys
    a{i} = permute(a{i}, [2,1,3]);
    a{i} = reshape(a{i}, n(i), r(i)*r(i+1));
    a{i} = fac{i}*a{i};
    a{i} = reshape(a{i}, size(fac{i},1), r(i), r(i+1));
    a{i} = permute(a{i}, [2,1,3]);
end
a = cell2core(tt_tensor, a);
a = full(a);
if ( nargin > 1 && ~isempty(sizes) )
    if (numel(sizes)==1)
        sizes = [sizes, 1];
    end;
   a=reshape(a,sizes);
end
return;
