function [a]=times(b,c,varargin)
%[A]=TIMES(B,C)
%[A]=TIMES(B,C,EPS)
%Hadamard product of two TT-tensors
%If optional parameter EPS is specified, use Krylov method
if (nargin == 2 )
    %Put zaglushka in it
%     a=tt_tensor(tt_hdm(core(b),core(c)));
%     return
    %Replace zaglushka by a code
    d=b.d; n=b.n; crb=b.core; crc=c.core; rb=b.r; rc=c.r;
    psb=b.ps; psc=c.ps;
    %Ranks are just a product
    a=tt_tensor;
    ra=rb.*rc;
    psa=cumsum([1;n(1:d).*ra(1:d).*ra(2:d+1)]);
    sz=dot(n(1:d).*ra(1:d),ra(2:d+1));
    cra=zeros(sz,1);
    for i=1:d
       cr1=crb(psb(i):psb(i+1)-1);
       cr2=crc(psc(i):psc(i+1)-1);
       cr1=reshape(cr1,[rb(i),n(i),rb(i+1)]);       
       cr2=reshape(cr2,[rc(i),n(i),rc(i+1)]);
       % rewrite in in the vectorized form
       cr1 = permute(cr1, [2, 1, 3]);
       cr1 = reshape(cr1, n(i)*rb(i)*rb(i+1), 1);
       ons = ones(1, rc(i)*rc(i+1));
       cr1 = cr1*ons;
       cr1 = reshape(cr1, n(i), rb(i), rb(i+1), rc(i), rc(i+1));
       cr1 = permute(cr1, [1, 2, 4, 3, 5]);
       cr1 = reshape(cr1, n(i), rb(i)*rc(i)*rb(i+1)*rc(i+1));
       cr2 = permute(cr2, [2, 1, 3]);
       cr2 = reshape(cr2, n(i)*rc(i)*rc(i+1), 1);
       ons = ones(1, rb(i)*rb(i+1));
       cr2 = cr2*ons;
       cr2 = reshape(cr2, n(i), rc(i), rc(i+1), rb(i), rb(i+1));
       cr2 = permute(cr2, [1, 4, 2, 5, 3]);
       cr2 = reshape(cr2, n(i), rb(i)*rc(i)*rb(i+1)*rc(i+1));       
       cr = cr1.*cr2;
       cr = reshape(cr, n(i), rb(i)*rc(i), rb(i+1)*rc(i+1));
       cr = permute(cr, [2, 1, 3]);
%        cr=zeros(ra(i),n(i),ra(i+1));
%        for j=1:n(i)
%           cr(:,j,:)=kron(reshape(cr1(:,j,:),rb(i),rb(i+1)),reshape(cr2(:,j,:),rc(i),rc(i+1)));
%        end
       cra(psa(i):psa(i+1)-1)=cr(:);
    end
    a.d=d;
    a.n=n;
    a.r=ra;
    a.ps=psa;
    a.core=cra;
elseif (nargin == 3)
     error('Krylov method is not yet supported, come again later');
end
return
end