function [p] = dot(qt1,qt2, do_qr)
%Dot product of two QTT-Tuckers
%   [PR]=DOT(TT1,TT2) -- dot product of two TT-tensors
%
%   [PR]=DOT(TT1,TT2, DO_QR) if DO_QR==true is specified, perform the 
%   left-to-right QRs of TT1,TT2
%   before the scalar product. It increases the  accuracy in some cases.
%
% In general, returns a 4D tensor of sizes 
% r0(tt1), r0(tt2), rd(tt1), rd(tt2)
% If r0(tt1) = r0(tt2) = 1 it returns a matrix of size rd(tt1) x rd(tt2)
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


if (nargin<3)||(isempty(do_qr))
    do_qr = false;
end;

d=qt1.dphys;
for i=1:d
    % dot product of factors gives rt1 times rt2 matrix, to be convolved
    % between cores
    % Make it more precise by QRs? <-- already inside @tt_tensor/dot
%     [qt1.tuck{i}, rv1]=qr(qt1.tuck{i}, 'lr');
%     [qt2.tuck{i}, rv2]=qr(qt2.tuck{i}, 'lr');
%     rv1 = eye(qt1.tuck{i}.r(end));
%     rv2 = eye(qt2.tuck{i}.r(end));
    Pfac = squeeze(dot(qt1.tuck{i}, qt2.tuck{i}, do_qr)); % size rt1,rt2
%     rt1 = size(rv1, 2); rt2 = size(rv2, 2);
%     rt1new = size(rv1,1); rt2new = size(rv2,1);
    % Now, merge Pfac to the core of qt1
    curcr = qt2.core{i};
    rc1 = size(curcr, 1); rt2 = size(curcr, 2); rc2 = size(curcr, 3);
    curcr = permute(curcr, [1, 3, 2]);
    curcr = reshape(curcr, rc1*rc2, rt2);
%     curcr = curcr*(rv1.');
    curcr = curcr*(Pfac.'); % Now, core2 has the tucker ranks of qt1
    curcr = reshape(curcr, rc1, rc2, qt1.core.n(i));
    qt2.core{i} = permute(curcr, [1, 3, 2]);
    
%     curcr = qt2.core{i};
%     rc1 = size(curcr, 1); rc2 = size(curcr, 3);
%     curcr = permute(curcr, [1, 3, 2]);
%     curcr = reshape(curcr, rc1*rc2, rt2);    
%     curcr = curcr*(rv2.');
%     curcr = reshape(curcr, rc1, rc2, rt2new);
%     qt2.core{i} = permute(curcr, [1, 3, 2]);    
end;
% Finaly, dot product of cores. It is consistent, since we merged Pfacs
p = dot(qt1.core, qt2.core, do_qr);
end  
