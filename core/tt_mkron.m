function [c]=tt_mkron(tt_arr)
%Compute Kronecker product of many TT-tensors. 
%   [C]=TT_MKRON(TT_ARR) Computes Kronecker product of many TT tensors
%   TT_ARR is a cell array of TT tensors, and the result RES is the 
%   Kronecker product of them. It is much faster, then taking them
%   one-by-one. 
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
N=numel(tt_arr);
df=zeros(N,1);
mm=zeros(N,1);
for i=1:N
  df(i)=ndims(tt_arr{i});
  mm(i)=mem(tt_arr{i});
end
full_mem=sum(mm);
full_d=sum(df);
cr=zeros(full_mem,1);
nf=zeros(full_d,1);
rf=zeros(full_d+1,1);
pos=1;
pos_d=1;
pos_r=1;
for i=1:N
    cur=tt_arr{i};
    cr(pos:pos+mm(i)-1)=cur.core(:);
    nf(pos_d:pos_d+df(i)-1)=cur.n(:);
    rf(pos_r:pos_r+df(i))=cur.r(:);
    pos=pos+mm(i);
    pos_d=pos_d+df(i);
    pos_r=pos_r+df(i);
end
c=tt_tensor;
c.core=cr;
c.n=nf;
c.r=rf;
c.d=full_d;
c.ps=cumsum([1;c.n.*c.r(1:c.d).*c.r(2:c.d+1)]);


return
end
