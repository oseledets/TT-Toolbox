function [tt]=round(tt,eps,rmax)
%[TT]=ROUND(TT,EPS)
%[TT]=ROUND(TT,EPS,RMAX)
%Approximate TT-matrix with relative accuracy EPS
if (nargin == 3 )
tt.tt=round(tt.tt,eps,rmax);
else
tt.tt=round(tt.tt,eps);    
end
return
end