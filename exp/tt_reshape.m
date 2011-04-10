function [res]=tt_reshape(tt,sz)
%[RES]=TT_RESHAPE(TT,SZ)
%Reshapes a TT-tensor into a new tensor.
%Only merging of dimensions is allowed right now!
d=sizel(tt,1);
if (size(sz,1) == 1)
  sz=sz';
end
dnew=size(sz,1);
res=cell(dnew,1);
%Decipher new sizes

return
end