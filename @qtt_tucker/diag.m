function [qt]=diag(qt)
% [QT]=DIAG(QT)
% Either makes a diagonal matrix from a vector in qtt_tucker, 
% or extracts a diagonal vector from a matrix

d = qt.dphys;
for i=1:d
    qt.tuck{i} = diag(qt.tuck{i});
end;

end