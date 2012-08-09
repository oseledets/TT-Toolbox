function ttf = friedman_test(n_problem, tt_rank, n_basis)
%FRIEDMAN regression

fun = @(i, j, x) x^(j-1);

if n_problem == 1
    %Friedman 1
    d = 10;
    r = [1; tt_rank*ones(d-1,1); 1];
    n = n_basis*ones(d,1);
    N = 20000;
    x = rand(d,N);
    y = zeros(1,N);
    for i = 1: N
        y(i) = 10*sin(pi*x(1,i)*x(2,i)) + 20*(x(3,i)-0.5)^2 + 10*x(4,i) +5*x(5,i);
    end
elseif n_problem == 2
    %Friedman 2
    d = 4;
    r = [1; tt_rank*ones(d-1,1); 1];
    n = n_basis*ones(d,1);
    N = 20000;
    x = rand(d,N);
    y = zeros(1,N);
    for i = 1: N
        x1 = x(1,i)*100;
        x2 = x(2,i)*520*pi + 40*pi;
        x3 = x(3,i);
        x4 = x(4,i)*10 + 1;
        y(i) = sqrt(x1^2 + (x2*x3 - 1/(x2*x4))^2);
    end
else
    %Friedman 3
    d = 4;
    r = [1; tt_rank*ones(d-1,1); 1];
    n = n_basis*ones(d,1);
    N = 20000;
    x = rand(d,N);
    y = zeros(1,N);
    for i = 1: N
        x1 = x(1,i)*100;
        x2 = x(2,i)*520*pi + 40*pi;
        x3 = x(3,i);
        x4 = x(4,i)*10 + 1;
        y(i) = atan((x2*x3 - 1/(x2*x4))/x1);
    end
end

var(y)
ttf = reg_als(x, y, fun, r, n)

end

