function [c]=kron2(a,b)
% Computes the Pseudokronecker product of TT tensors:
%   [C]=KRON2(A,B)
%   c{i}=kron(a{i},b{i}) <- "a" indices vary faster!
%   mode sizes and ranks are multiplied
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

if ( isempty(a) )
  c=b;
  return
elseif ( isempty(b) )
  c=a;
  return
end
if ( size(a.core,1) == 1 ) 
   a.core=(a.core).';
end
if ( size(b.core,1) == 1 ) 
   b.core=(b.core).';
end

d = a.d;
if (b.d~=d)
    error('Could not compute kron2 from the tensors with not equal dimensions\n');
end;
cra = a.core;
na = a.n;
ra = a.r;
psa = a.ps;
crb = b.core;
nb = b.n;
rb = b.r;
psb = b.ps;

c=tt_tensor;
n = [(na).*(nb)];
r = [(ra).*(rb)];
ps = cumsum([1; n.*r(2:d+1).*r(1:d)]);


c.d=d;
c.n=n;
c.r=r;
c.ps = ps;
cr = zeros(ps(d+1)-1, 1);
for i=1:d
    a1 = cra(psa(i):psa(i+1)-1);
    b1 = crb(psb(i):psb(i+1)-1);
    c1 = a1*(b1.'); % Size ra1*na*ra2, rb1*nb*rb2
    c1 = reshape(c1, ra(i), na(i), ra(i+1), rb(i), nb(i), rb(i+1));
    c1 = permute(c1, [1, 4, 2, 5, 3, 6]); % We need rc1, nc*rc2
    cr(ps(i):ps(i+1)-1) = c1(:);
end;

c.core=cr;

end