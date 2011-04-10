%Converts a linear index to a multiindex
function [sub]=ind_to_sub(ind,sz)
if ( size(sz,2) == 1 )
   d=size(sz,1);
   sz=sz';
else
    d=size(sz,2);
end
sub=zeros(1,d);
sub=ind2sub(sz,ind);
%d=size(sz,2);
%ind = i1+(j1-1)*sz(2)+(k1-1)*sz(2)*sz(3) + ... + 
%prs=1;
%mult=1;
%for k=1:d
%    mult=mult+(ind-1)*prs;
%    prs = prs*sz(k);
%end
return