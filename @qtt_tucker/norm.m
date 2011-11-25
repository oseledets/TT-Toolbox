function [nrm] = norm(tt)
%[NRM]=NORM(QTT_TUCKER)
%Frobenius norm of the QTT_TUCKER
d=tt.dphys;
core=tt.core;
tuck=tt.tuck;
for i=1:d
   [tuck{i},rm]=qr(tuck{i},'lr');
   core{i}=ten_conv(core{i},2,rm.');
end
nrm=norm(core);

return
end