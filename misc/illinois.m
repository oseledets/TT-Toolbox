function x=illinois(f,a,b,tol)
fa=f(a); fb=f(b);
while abs(b-a)>tol
    step=fa*(a-b)/(fb-fa);
    x=a+step;
    fx=f(x);
    if sign(fx)~=sign(fb)
        a=b; fa=fb;
    else
        fa=fa/2;
    end
    b=x; fb=fx;
end