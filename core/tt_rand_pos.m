function [tt]=tt_rand_pos(n, d, r)
% function [tt]=tt_rand_pos(n, d, r)
% Generates POSITIVE (important for quadr_solve) pseudo-random tensor

if (max(size(n))==1)
    n = n*ones(d,1);
end;
if (max(size(r))==1)
    r = r*ones(d-1,1);
end;

r = [1; r; 1];
ps = cumsum([1; r(1:d).*n.*r(2:d+1)]);
cr = zeros(ps(d+1)-1,1);

for i=1:d
    curcore = rand(r(i), n(i), r(i+1));
    cr(ps(i):ps(i+1)-1)=curcore(:);
end;

tt = tt_tensor;
tt.d = d;
tt.n = n;
tt.r = r;
tt.ps = ps;
tt.core = cr;

end