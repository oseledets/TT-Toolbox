function [tt1]=conj(tt)
%[TT1]=CONJ(TT)
%Compute complex conjugate of TT-tensor TT
tt1=tt;
tt1.core=conj(tt.core);
return
end