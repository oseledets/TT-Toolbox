function a=tt_dum1(a,b)
%Technical function. Never try to know what (and why!) it does
d=a.d; r=a.r; n=a.n;
a{a.d}=reshape(reshape(a{d},[r(d)*n(d),r(d+1)])*b,[r(d),n(d),size(b,2)]);
return
end
