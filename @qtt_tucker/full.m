function [a] = full(tt, sizes)
%[A]=FULL(QTT_TUCKER)
%Converts a QTT-Tucker tensor the full tensor
a=full(tt.core); 
fac=cell(tt.dphys,1);
for i=1:tt.dphys
   ti=tt.tuck{i}; ri=rank(ti,ti.d+1);
    ti=full(ti);
   fac{i}=reshape(ti,[numel(ti)/ri,ri]);
end
%And now convert tucker representation to the full one
n=size(tt.core);
for i=1:tt.dphys
   a=reshape(a,[n(i),numel(a)/n(i)]);
   a=(fac{i}*a).'; %r1r2r3 -> n1r2r3 ->r2r3n1 
end
if ( nargin > 1 && ~isempty(sizes) )
   a=reshape(a,sizes);
end
return;
