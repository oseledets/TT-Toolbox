function [res]=tt_smv(sttm, vec)
%TT-matrix with sparse factors by TT-vector multiplication
%   [res]=TT_SMV(STTM,VEC) Multiplies the TT-matrix stored in the TT1.0
%   format but with "sparse" cores by a vector VEC
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


d=size(vec,1);
res=cell(d,1);

n=sttm{d+2}(1);
% m=size(vec{1},1);
rm1=sttm{d+1}(1);
rv1=size(vec{1},2);
%mat{1} is (n x m x rm1 , vec{1} is mxrv1)
%            n x rm1 x rv1
res{1}=sttm{1}*vec{1};
res{1}=reshape(permute(reshape(res{1},[n,rm1,rv1]),[1,3,2]),[n,rv1*rm1]);

if (d>1)
    n=sttm{d+2}(d);
    % m=size(mat{d},2);
    rm1=sttm{d+1}(d-1);
    rv1=size(vec{d},2);
    %mat{1} is (n x rm1 x rv1 , vec{1} is mxrv1)
    res{d}=sttm{d}*vec{d};
    res{d}=reshape(permute(reshape(res{d},[n,rm1,rv1]),[1,3,2]),[n,rv1*rm1]);
end;

for i=2:d-1
    n=sttm{d+2}(i);
    m=size(vec{i},1);
    rm1=sttm{d+1}(i-1);
    rm2=sttm{d+1}(i);
    rv1=size(vec{i},2);
    rv2=size(vec{i},3);
    %mat{i} is n x m x rm1 x rm2, vec{i} is m x rv1 x rv2
    %n x (rv1*rv2)*rm1*rm2
    %n x rm2 x rm1 x m, n x rm2 x rm1 x (rv1*rv2)
    %n x rm2 x rmx1 x rv1 x rv2
    %want: n x rv1*rm1 x rv2*rm2
    res{i}=sttm{i}*reshape(vec{i},[m,rv1*rv2]);
    res{i}=reshape(res{i},[n,rm2,rm1,rv1,rv2]);
    res{i}=permute(res{i},[1,4,3,5,2]);
    res{i}=reshape(res{i},[n,rv1*rm1,rv2*rm2]);
    %res{i}=permute(reshape(permute(mat{i},[1,3,4,2])*reshape(vec{i},[m,rv1*rv2]),[n,rv1,rv2,rm1,rm2]),[1,2,4,3,5]);
end 

end