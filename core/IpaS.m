function [M]=IpaS(d, a)
% function [M]=IpaS(d, a)
% Generates I+a*S_{-1} matrix in the QTT format:
%  1 0 0 0
%  a 1 0 0
%  0 a 1 0
%  0 0 a 1
% Convenient for Crank-Nicolson and time gradient matrices

M = cell(d,1);
M{1} = zeros(2,2,2);
M{1}(:,:,1)=[1,0;a,1];
M{1}(:,:,2)=[0,a;0,0];
for i=2:d-1
    M{i}=zeros(2,2,2,2);
    M{i}(:,:,1,1)=eye(2);
    M{i}(:,:,2,1)=[0,0;1,0];
    M{i}(:,:,2,2)=[0,1;0,0];
end;
M{d} = zeros(2,2,2);
M{d}(:,:,1)=eye(2);
M{d}(:,:,2)=[0,0;1,0];

end