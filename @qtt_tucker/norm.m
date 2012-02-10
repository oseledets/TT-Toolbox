function [nrm] = norm(tt)
%Frobenius norm of the QTT-Tucker
%   [NRM]=NORM(QTT_TUCKER) Computes the Frobenius norm of the QTT-Tucker
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
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