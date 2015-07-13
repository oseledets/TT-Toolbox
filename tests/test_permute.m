rng(42);
for i = 1:1000
    d = randi(9)+1;
    n = randi(4,[1,d]);
    x = tt_rand(n,d,[1,randi(20,[1,d-1]),1]);
    order = randperm(d);
    tol = 10^(-5*rand());

    yf = permute(reshape(full(x), n),order);
    yn = reshape(full(permute(x, order, tol)),n(order));
    err = yf - yn;
    err = norm(err(:))/norm(yf(:));
    if (err > tol) warning(sprintf('Error (%g) is larger than tolerance (%g)', err, tol)); end;
end

