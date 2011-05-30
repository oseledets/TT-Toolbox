d = 12;

a = 0;
b = 2;

h = (b-a)/(2^d+1);

x = (1:1:2^d)'*h;
x = tt_tensor(full_to_qtt(x,1e-12));

e = tt_tensor(tt_ones(d,2));

y = kron(e,x);
x = kron(x,e);
x2 = x.*x;
x2 = round(x2, 1e-12);

f = (kron(e,e)-x).*(kron(e,e)-x);
f = f+100*(y-x2).*(y-x2);
f = round(f, 1e-12);

S1 = dct(eye(16));
I1 = tt_matrix(tt_eye(2,1));
S1 = tt_matrix(full_to_nd(S1, 1e-12, 4));
S = tt_matrix(tt_eye(2,2*d));
for i=1:2:2*d-3
    cur_S = S1;
    for j=1:i-1
        cur_S = kron(I1, cur_S);
    end;
    for j=i+1:2*d-3
        cur_S = kron(cur_S, I1);
    end;
    S = S*cur_S;
end;


F = diag(f);

F2 = round(S*F, 1e-6)*S';
F2 = round(F2, 1e-6);

u = dmrg_eig2(core(F2), tt_random(2,2*d,2), 1e-6, 150);
u = tt_tensor(u)

