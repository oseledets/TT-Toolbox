function [p] = dot(qt1,qt2)
%[P]=DOT(QT1,QT2)
%Dot  product of two QTT-tuckers
%
%

d=qt1.dphys;
for i=1:d
    % dot product of factors gives rt1 times rt2 matrix, to be convolved
    % between cores
    % Make it more precise by QRs?
    [qt1.tuck{i}, rv1]=qr(qt1.tuck{i}, 'lr');
    [qt2.tuck{i}, rv2]=qr(qt2.tuck{i}, 'lr');
%     rv1 = eye(qt1.tuck{i}.r(end));
%     rv2 = eye(qt2.tuck{i}.r(end));
    Pfac = dot(qt1.tuck{i}, qt2.tuck{i});
    rt1 = size(rv1, 2); rt2 = size(rv2, 2);
    rt1new = size(rv1,1); rt2new = size(rv2,1);
    % Now, merge Pfac to the core of qt1
    curcr = qt1.core{i};
    rc1 = size(curcr, 1); rc2 = size(curcr, 3);
    curcr = permute(curcr, [1, 3, 2]);
    curcr = reshape(curcr, rc1*rc2, rt1);
    curcr = curcr*(rv1.');
    curcr = curcr*Pfac; % Now, core1 has the tucker ranks of qt2
    curcr = reshape(curcr, rc1, rc2, rt2new);
    qt1.core{i} = permute(curcr, [1, 3, 2]);
    
    curcr = qt2.core{i};
    rc1 = size(curcr, 1); rc2 = size(curcr, 3);
    curcr = permute(curcr, [1, 3, 2]);
    curcr = reshape(curcr, rc1*rc2, rt2);    
    curcr = curcr*(rv2.');
    curcr = reshape(curcr, rc1, rc2, rt2new);
    qt2.core{i} = permute(curcr, [1, 3, 2]);    
end;
% Finaly, dot product of cores. It is consistent, since we merged Pfacs
p = dot(qt1.core, qt2.core);
end  
