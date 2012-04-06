function [p] = dot(tt1,tt2,do_qr)
%Dot  product of two TT tensors
%   [PR]=DOT(TT1,TT2) -- dot product of two TT-tensors
%
%   [PR]=DOT(TT1,TT2, DO_QR) if DO_QR==true is specified, perform the 
%   left-to-right QRs of TT1,TT2
%   before the scalar product. It increases the  accuracy in some cases.
%
% In general, returns a 4D tensor of sizes 
% r0(tt1), r0(tt2), rd(tt1), rd(tt2)
% If r0(tt1) = r0(tt2) = 1 it returns a matrix of size rd(tt1) x rd(tt2)
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


if (nargin<3)||(isempty(do_qr))
    do_qr = false;
end;

if (do_qr)
    [tt1,rv1]=qr(tt1, 'lr');
    [tt2,rv2]=qr(tt2, 'lr');
end;

d=tt1.d; 
r1=tt1.r; r2=tt2.r; ps1=tt1.ps; ps2=tt2.ps;
n=tt1.n;
core1=tt1.core; core2=tt2.core;

%ps is always r1(i-1)xr; but if there is a hanging thing? what to do?
%That means, I define a dot product of two "hanging" tensors as a matrix...
%brr.... And if it is hanging on the right? 
% 
% p=ones(r1(1),r2(1)); % Over r1(1) and r2(1) there is not summation blin.
% %So, if we just sum over n(1) separatedly and leave a strange QxR1(I)xR2(I)
% %matrix... 
% for i=1:d
%   cr1=core1(ps1(i):ps1(i+1)-1);
%   cr2=core2(ps2(i):ps2(i+1)-1);
%   p=reshape(p,[r1(i),r2(i)]);
%   cr2=reshape(cr2,[r2(i),numel(cr2)/r2(i)]);
%   p=p*cr2; %p is Q*r1(i)xn(i)xr2(i+1);
%   cr1=reshape(cr1,[r1(i)*n(i),numel(cr1)/(r1(i)*n(i))]);
%   p=reshape(p,[r1(i)*n(i),numel(p)/(r1(i)*n(i))]);
%   p=cr1'*p;
% end

% If the first indices are not ones
p=eye(r1(1)*r2(1));
p = reshape(p, r1(1)*r2(1)*r1(1), r2(1));

for i=1:d
    cr1=core1(ps1(i):ps1(i+1)-1);
    cr2=core2(ps2(i):ps2(i+1)-1);
    cr2=reshape(cr2,[r2(i), n(i)*r2(i+1)]);
    
    p = p*cr2; % size r11*r21*r1-, n*r2+
    p = reshape(p,r1(1)*r2(1), r1(i)*n(i), r2(i+1));
    p = permute(p, [1, 3, 2]);
    p = reshape(p, r1(1)*r2(1)*r2(i+1), r1(i)*n(i));
    
    cr1=reshape(cr1,[r1(i)*n(i), r1(i+1)]);
    
    p = p*conj(cr1); % size r11*r12*r2+, r1+
    p = reshape(p, r1(1)*r2(1), r2(i+1), r1(i+1));
    p = permute(p, [1, 3, 2]);
    p = reshape(p, r1(1)*r2(1)*r1(i+1), r2(i+1));
end;

if (do_qr)
    r2old = size(rv2, 2);
    r1old = size(rv1,2);
    p = p*rv2;
    p = reshape(p, r1(1)*r2(1), r1(d+1), r2old);
    p = permute(p, [1, 3, 2]);
    p = reshape(p, r1(1)*r2(1)*r2old, r1(d+1));
    p = p*conj(rv1);
    p = reshape(p, r1(1), r2(1), r2old, r1old);
    p = permute(p, [1,2,4,3]);
    if ( r1(1)*r2(1) == 1 ) %Save the rabbits
        p=reshape(p, r1old, r2old);
    end
else
    p = reshape(p, r1(1), r2(1), r1(d+1), r2(d+1));
    if ( r1(1)*r2(1) == 1 ) %Save the rabbits
        p=reshape(p, r1(d+1), r2(d+1));
    end
end;
end
