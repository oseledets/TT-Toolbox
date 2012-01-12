function [tt]=ctranspose(tt)
%[TT1]=CTRANSPOSE(TT)
%Compute complex conjugate transpose of a TT-matrices insite qtt_tucker

d = tt.dphys;
for i=1:d
    tt.tuck{i} = ctranspose(tt.tuck{i});
end;