function test_tt_dot_complex()
% From the toolbox root: setup; addpath('tests'); test_tt_dot_complex
e = [1; 0];
x = {e; e};
y = {1i*e; -1i*e};
% Both tensor products represent [1;0;0;0]. The old contraction returns -1.
assert(abs(tt_dot(x,y)-1) < 1e-14);
assert(abs(tt_dot({[1;1i]},{[1;1]})-(1-1i)) < 1e-14);

% Exercise an interior core and compare to the dense Hermitian product.
j = (1:8).';
a = sin(j)+1i*cos(2*j);
b = cos(j.^2)-1i*sin(j/3);
u = core(tt_tensor(reshape(a,[2 2 2]),1e-14));
v = core(tt_tensor(reshape(b,[2 2 2]),1e-14));
assert(abs(tt_dot(u,v)-a'*b) < 1e-12*norm(a)*norm(b));
end
