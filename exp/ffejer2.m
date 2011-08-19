function [x,w]=ffejer2(N, a, b)

N2 = N+1;
x = a+(b-a)*0.5*(cos(pi*(0:N2)/N2)+1)';
x = x(2:N+1);
w = zeros(N2,1);

for i=0:N2
    curw = 0;
    for j=1:floor(N2/2)
        curw = curw + sin((2*j-1)*i*pi/N2)/(2*j-1);
    end;
    w(i+1)=curw*4/N2*sin(i*pi/N2);
end;
w = w(2:N+1);
w = w*(b-a)/2;

end