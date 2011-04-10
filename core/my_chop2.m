 function [r] = my_chop2(sv,eps)
%[R]=MY_CHOP2(SV,EPS) 
%Truncation by absolute precision in Frobenius norm
%Find minimal possible R such that \sqrt(sum(SV(R:n)^2)) < EPS
 r = 0;
 if ( size(sv,1) == 1 )
   sv=sv';
 end
 n=size(sv,1);
 er=norm(sv);
 while ( (er > eps) && ( r < n ))
   r=r+1;
   er=norm(sv((r+1):n));
   %er = er - sv(r)^2;
 end
 if ( r == (n+1) )
   r=n;
 end   
 return
 end
