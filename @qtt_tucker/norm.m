function [nrm] = norm(tt)
%[NRM]=NORM(QTT_TUCKER)
%Frobenius norm of the QTT_TUCKER
d=tt.dphys;
core=tt.core;
tuck=tt.tuck;
for i=1:d
    if (isa(tuck{i}, 'tt_matrix'))
        tuck{i} = tt_tensor(tuck{i});
    end;
   [tuck{i},rm]=qr(tuck{i},'lr');
   core{i}=ten_conv(core{i},2,rm.');
end
nrm=norm(core);

return
end