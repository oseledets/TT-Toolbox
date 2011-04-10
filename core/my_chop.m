function [r]=my_chop(a,eps)
%[R]=MY_CHOP(A,EPS)
%Chops values of a vector A by relative threshold EPS
nrm2=norm(a).^2;
bound=eps.^2*nrm2;
[r1,tmp]=size(a(:));
er=0;
r=r1;
while(er < bound && r >= 1 )
    er=er+a(r).^2;
    r=r-1;
end
r=r+1;
return
end
