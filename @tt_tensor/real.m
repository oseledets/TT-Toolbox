function [b]=real(a)
%Compute a real part of TT-tensor 
%   [B]=REAL(A)
%
% TT-Toolbox 2.2, 2009-2013
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%
% For details see Thm4 in Dolgov, Khoromskij, Savostyanov,
% "Superfast Fourier transform using the QTT format"
% http://dx.doi.org/10.1007/s00041-012-9227-4
%
%---------------------------

if (a.d==1)
    b=a;
    b.core=real(b.core);
    exit;
end
d=a.d;

%Determine the size and preallocate the result
b=tt_tensor;
b.r=a.r; b.r(2:d)=2*a.r(2:d);
b.d=a.d; b.n=a.n;
sz=dot(b.n.*b.r(1:d),b.r(2:d+1));
pos=(b.n.*b.r(1:d)).*b.r(2:d+1);
b.ps=cumsum([1;pos]);
b.core=zeros(sz,1);
apos=a.ps; bpos=b.ps;

% Fill cores
for i=1:d
    bb=zeros(b.r(i),b.n(i),b.r(i+1));
    aa=reshape(a.core(apos(i):apos(i+1)-1), [a.r(i),a.n(i),a.r(i+1)]);
    if (i==1)
        bb(1:a.r(i),:,1:a.r(i+1)           )=real(aa);
        bb(1:a.r(i),:,1+a.r(i+1):2*a.r(i+1))=imag(aa);
    elseif (i==d)
        bb(1:a.r(i),         :,1:a.r(i+1))= real(aa);
        bb(1+a.r(i):2*a.r(i),:,1:a.r(i+1))=-imag(aa);
    else
        bb(1:a.r(i),         :,1:a.r(i+1)           )= real(aa);
        bb(1:a.r(i),         :,1+a.r(i+1):2*a.r(i+1))= imag(aa);
        bb(1+a.r(i):2*a.r(i),:,1:a.r(i+1)           )=-imag(aa);
        bb(1+a.r(i):2*a.r(i),:,1+a.r(i+1):2*a.r(i+1))= real(aa);
    end
    b.core(bpos(i):bpos(i+1)-1)=bb(:);
end
return
end