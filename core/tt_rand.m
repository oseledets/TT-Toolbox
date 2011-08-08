function [tt]=tt_rand(n,d,r)
%[TT]=TT_RAND(N,D,R)
%Generate a random tensor
%N --- is a mode size array
%It can be either a number, or an array of dimension d
%R --- rank array.
%Either a number, or an array of dimension d+1

if ( numel(n) == 1 ) 
 n = n * ones(d,1);
end
if ( numel(r) == 1 )
 r = r * ones(d+1,1); 
 r(1) = 1;
 r(d+1)=1;
end
tt=tt_tensor;
tt.n=n;
tt.r=r;
tt.d=d;
ps=cumsum([1;n.*r(1:d).*r(2:d+1)]);
tt.ps=ps;
cr=randn(ps(d+1)-1,1);
tt.core=cr;
return
end