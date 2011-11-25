function [c]=tt_dum1(a,b)
    c=a;
    d=a.d;
    n=a.n;
    r=a.r;
    mm=a{d};
    mm=reshape(mm,[r(d)*n(d),r(d+1)]);
    mm=mm*b;
    c{d}=reshape(mm, r(d), n(d), size(b,2));
return
end