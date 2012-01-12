function [tt]=tt_qreshape(tt,s,sz)

% Reshapes the cores of the input (Q)TT representation.
% s (scalar) is the number of mode indices in each core of the input representation;
% sz is a matrix of size d x s1, where d is the number of cores and s1 is the number of mode indices in each core of the output representation;
% sz(k,alpha) is the mode size of the alpha-th index in the k-th core of the output representation.
% k-th component of prod(sz,2) is the overall numbers of mode degrees of freedom of the k-th core of the output representation 
%
% April 26, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

d=size(tt,1);
for k=1:d
	szk=size(tt{k});
	szk=[szk,ones(1,s)];
	szkr=szk(s+1:numel(szk));
	tt{k}=reshape(tt{k},[sz(k,:),szkr]);
end

return
end