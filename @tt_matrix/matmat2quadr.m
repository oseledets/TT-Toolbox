function [Q]=matmat2quadr(A, B)
% function [Q]=matmat2quadr(A, B)
% Computes 3D tt_array from tt_matrices:
% Q_{i,j,k} = A_{i,j} B_{i,k}

ttA = A.tt;
na = A.n;
ma = A.m;
crA = ttA.core;
rA = ttA.r;
psA = ttA.ps;
d = ttA.d;

ttB = B.tt;
nb = B.n;
mb = B.m;
crB = ttB.core;
psB = ttB.ps;
rB = ttB.r;

if (na~=nb)
    fprintf('The first mode sizes shall be equal! Consider using (c)transpose\n');
    return;
end;

Q = tt_array;
Q.ns = [na,ma,mb];
% Q.m = ma;
% Q.k = mb;

rQ = rA.*rB;
psQ = zeros(d+1,1);
% Preallocate Q-core
crQ = zeros(sum(rQ(1:d).*na.*ma.*mb.*rQ(2:d+1)), 1);
psQ(1)=1;

for i=1:d
    a1 = crA(psA(i):psA(i+1)-1);
    a1 = reshape(a1, rA(i), na(i), ma(i), rA(i+1));
    b1 = crB(psB(i):psB(i+1)-1);
    b1 = reshape(b1, rB(i), nb(i), mb(i), rB(i+1));
    
    q1 = zeros(rA(i)*rB(i), na(i), ma(i), mb(i), rA(i+1)*rB(i+1));
    for j=1:na(i)
        cur_a1 = reshape(a1(:,j,:,:), rA(i), ma(i)*rA(i+1));
        cur_b1 = reshape(b1(:,j,:,:), rB(i), mb(i)*rB(i+1));
        cur_q1 = kron(cur_b1, cur_a1); % size rA1*rB1, ma1*rA1*mb1*rB1
        cur_q1 = reshape(cur_q1, rA(i)*rB(i), ma(i), rA(i+1), mb(i), rB(i+1));
        cur_q1 = permute(cur_q1, [1 2 4 3 5]);
        cur_q1 = reshape(cur_q1, rA(i)*rB(i), 1, ma(i), mb(i), rA(i+1)*rB(i+1));
        q1(:,j,:,:,:)=cur_q1;
    end;
    
    crsize = rQ(i)*na(i)*ma(i)*mb(i)*rQ(i+1);
    crQ(psQ(i):(psQ(i)+crsize-1)) = reshape(q1, crsize, 1);
    psQ(i+1) = psQ(i)+crsize;
end;

ttQ = tt_tensor;
ttQ.d = d;
ttQ.core = crQ;
ttQ.n = na.*ma.*mb;
ttQ.r = rQ;
ttQ.ps = psQ;

Q.tt = ttQ;

end