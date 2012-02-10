function [res]=tt_dist3(tt1,tt2)
%Accurate distance between two tensors in TT1.0 format
%   [RES]=TT_DIST3(TT1,TT2) More accurate and but more expensive than
%   TT_DIST computation of a distance between two tensors
%   TT1 and TT2 in Frobenius norm: just recompression
%   of their difference. Please avoid its usage: it will be removed in
%   future releases. Use norm(a-b) from the object-oriented version
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
rs=tt_add(tt1,tt_scal(tt2,-1.0));
rs{1}=reshape(rs{1}, size(rs{1},1),1,size(rs{1},2));
d=size(rs,1);
for i=d:-1:2
    cr = rs{i}; n1 = size(cr,1); r1 = size(cr,2); r2 = size(cr,3);
    cr = reshape(permute(cr, [1 3 2]), n1*r2, r1);
    [cr,rv]=qr(cr,0);
    nrm1 = norm(rv, 'fro');
    cr = cr*nrm1; rv = rv/nrm1;
    r = size(cr,2);
    rs{i}=permute(reshape(cr, n1, r2, r), [1 3 2]);
    cr2 = rs{i-1}; n0 = size(cr2,1); r0 = size(cr2, 2);
    cr2 = reshape(cr2, n0*r0, r1);
    cr2 = cr2*(rv.');
    rs{i-1}=reshape(cr2, n0, r0, r);
end;
rs{1}=reshape(rs{1}, size(rs{1},1),size(rs{1},3));
% rs=tt_compr2(rs,0.d0); %After recompression cores are - Are you crazy,
% dude?!?! Who will make a full SVD just to compute a norm?!
res=norm(rs{1}(:))/sqrt(size(rs{1},2));
res=res*norm(rs{d}(:));
for k=2:d-1
  res = res * norm(rs{k}(:))/sqrt(size(rs{k},3));
end
return
end
