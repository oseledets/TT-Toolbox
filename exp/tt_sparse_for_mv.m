function [sttm] = tt_sparse_for_mv(ttm, cuttol)
%Sparsify the TT matrix in TT1.0 matrix format
%   [STTM]=TT_SPARSE_FOR_MV(TTM,CUTTOL) Computes the sparse version of the
%   matrix in the TT1.0 matrix format. May work for certain operators
%
%
% TT Toolbox 2.1, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
z
d=size(ttm,1);

sttm=cell(d+2,1);
sttm{d+1}=zeros(d-1,1);
sttm{d+2}=zeros(d,1);

n=size(ttm{1},1);
m=size(ttm{1},2);
rm1=size(ttm{1},3);
sttm{d+1}(1)=rm1;
sttm{d+2}(1)=n;
sttm{1}=reshape(permute(ttm{1},[1,3,2]),[n*rm1,m]);
sttm{1}(find(abs(sttm{1})/max(max(abs(sttm{1})))<cuttol))=0;
sttm{1}=sparse(sttm{1});

n=size(ttm{d},1);
m=size(ttm{d},2);
rm1=size(ttm{d},3);
sttm{d+1}(d-1)=rm1;
sttm{d+2}(d)=n;
sttm{d}=reshape(permute(ttm{d},[1,3,2]),[n*rm1,m]);
sttm{d}(find(abs(sttm{d})/max(max(abs(sttm{d})))<cuttol))=0;
sttm{d}=sparse(sttm{d});

for i=2:d-1
    n=size(ttm{i},1);
    m=size(ttm{i},2);
    rm1=size(ttm{i},3);
    rm2=size(ttm{i},4);
    sttm{d+1}(i)=rm2;    
    sttm{d+2}(i)=n;
    sttm{i}=reshape(permute(ttm{i},[1,4,3,2]),[n*rm2*rm1,m]);
    sttm{i}(find(abs(sttm{i})/max(max(abs(sttm{i})))<cuttol))=0;
    sttm{i}=sparse(sttm{i});
end 

end