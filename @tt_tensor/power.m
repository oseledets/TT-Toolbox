function [c] = power(a,b,varargin)
%[C]=POWER(A,B)
%[C]=POWER(A,B,EPS)
%Compute A.^B for a TT-tensor A and natural B.
%Optional accuracy parameter EPS --- can use Krylov method in future 
if isa(a,'tt_tensor') && numel(b) == 1
    c=a;
    b=b-1;
    while(b>0)
      c=c.*a;
      b=b-1;
    end
    
else
    error('Incorrect usage for TT-power');
end
