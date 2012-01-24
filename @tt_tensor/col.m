function [y]=col(x,k)
%[Y]=COL(X,K)
%Computes the subcolumns of a 
%It is useful for block algorithms when r(d+1) is not equal to 1
%For analogous thing on r(1) see ROW(X,K)

%d=x.d;
%r=x.r;
%n=x.n;
%cr=x.core;
%ps=x.ps;
%y=tt_tensor;

d=x.d;
rx=x.r;
nx=x.n;
cr=subsref(x,substruct('{}',{d}));
cr=reshape(cr,[numel(cr)/rx(d+1),rx(d+1)]);
cr=cr(:,k);
cr=reshape(cr,[rx(d),nx(d),numel(k)]);
%x{d}=reshape(x{d},[rx(d),nx(d),numel(k)]);
x=subsasgn(x,substruct('{}',{d}),cr);
y=x;

%Modification is required only for the last core
%cr_last=cr(ps(d):ps(d+1)-1);
%cr_last=reshape(cr_last,[r(d),n(d),r(d+1)]);
%cr_last=cr_last(:,:,k);
%r(d+1)=size(cr_last,3);
%cr(ps(d):ps(d+1)-1)=[];
%cr=[cr,cr_last(:)'];
%r(d+1)=size(cr_last,3);
%ps(d+1)=ps(d)+r(d)*n(d)*r(d+1);
%y.n=n;
%y.d=d;
%y.r=r;
%y.ps=ps;
%y.core=cr;
%return
end
