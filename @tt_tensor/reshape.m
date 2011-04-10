function [tt1]=reshape(tt,sz,varargin)
%[TT1]=RESHAPE(TT,SZ)
%[TT1]=RESHAPE(TT,SZ,EPS)
%Reshapes TT-tensor into a new one, with dimensions specified by SZ
%Optionally, accuracy EPS can be specified
d=tt.d;
n=tt.n;
ps=tt.ps;
cr=tt.core;
d0=numel(sz);
if ( prod(sz) ~= prod(n) )
 error('Reshape: incorrect sizes \n');
end
% tt1=tt_tensor;
% %Now we have to determine the set of elementary operations.
% %In new tensor, dimensions are either merged, or split, and this
% %'diagram' has to be derived
% i1=1;
% cr1=[];
% 
% for i=1:d
%   %Determine what to do with the first core: split or merge
%   if ( sz(i1) < n(i) ) %split
%       %Determine the number of cores to split
%       cm=cumprod(sz); 
%       ff=find(cm>n(i)); ff=ff(1); ff=ff-1;
%       spt=sz(1:ff); %These are new! dimensions
%       core=cr(ps(i):ps(i+1)-1); %Core to split
%       core=reshape(core,[r(i),spt
%   elseif ( sz(i1) == n(i) ) %do nothing
%   else %Merge
%       %Determine the number of cores to merge
%   end
% end
% if ( nargin == 2 )
%     
% elseif ( nargin == 3 )
%     
% end
    error('Reshape function is not implemented yet');
return
end