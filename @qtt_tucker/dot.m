function [p] = dot(qt1,qt2)
%[P]=DOT(QT1,QT2)
%Dot  product of two QTT-tuckers
%
%

d=qt1.dphys;
for i=1:d
    % dot product of factors gives rt1 times rt2 matrix, to be convolved
    % between cores
    Pfac = dot(qt1.tuck{i}, qt2.tuck{i});
    rt1 = size(Pfac, 1); rt2 = size(Pfac, 2);
    % Now, merge Pfac to the core of qt1
    curcr = qt1.core{i};
    rc1 = size(curcr, 1); rc2 = size(curcr, 3);
    curcr = permute(curcr, [1, 3, 2]);
    curcr = reshape(curcr, rc1*rc2, rt1);
    curcr = curcr*Pfac; % Now, core1 has the tucker ranks of qt2
    curcr = reshape(curcr, rc1, rc2, rt2);
    qt1.core{i} = permute(curcr, [1, 3, 2]);
end;
% Finaly, dot product of cores. It is consistent, since we merged Pfacs
p = dot(qt1.core, qt2.core);
end  
