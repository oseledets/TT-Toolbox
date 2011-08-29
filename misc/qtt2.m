%Test a qtt2 idea by Oswald

a=0;
b=1;
d=10;
n=2^d-1;
h=(b-a)/(n-1);
x=a+(0:n-1)*h;
%Create a sine function with noise

f=sin(10*x); 
f=f+0*1e-1*randn(size(f));
plot(f)

%Now do QTT2

f1=[f(1:2:n-1);f(2:2:n-1);f(3:2:n)];

%1 2 3  
%2 3 4  
%3 4 5 

%1 5   9            1 3 5 5 7  9   9  11 13
%2 6   10           2 4 6 6 8  10  10 12 14
%3 7   11           3 5 7 7 9  11  11 13 15
%4 8   12  ->       4 6 8 8 10 12  12 14 16
%5 9   13           
%6 10  14
%7 11  15
%8 12  16


% 1 3 5
% 2 4 6
% 3 5 7
% 4 6 8






