function [c] = mtimes(a,b,varargin)

if (isa(a,'tt_array') && isa(b,'double') && numel(b) == 1 || (isa(b,'tt_array') && isa(a,'double') && numel(a) == 1) )
%Scalar    
   c=tt_quadr;
    if (  isa(b,'double') )
      c.ns=a.ns; c.tt=(a.tt)*b;
   else
      c.ns=b.ns; c.tt=(b.tt)*a;
    end
end;

   
end
