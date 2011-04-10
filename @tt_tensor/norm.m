function [nrm] = norm(tt)
%[NRM]=NORM(TT)
%Frobenius norm of the tt_tensor
%Zatycka

nrm=sqrt(abs(dot(tt,tt)));
return
end