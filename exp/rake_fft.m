function [y]=rake_fft(x, tol)
% Discrete Fourier Transform in the QTT-Tucker format
%   function [y]=rake_fft(x, tol)
%   Computes the multidimensional DFT of the vector x in the QTT-Tucker format with
%   the approximation tolerance tol.
%   Returns the quantized dimensions in the correct order.
%
%   !!WARNING!! employs a (probably) costly rank reverse procedure.

D = size(x.tuck);
tuck = x.tuck;
crx = core2cell(x.core);
for i=1:D
    [tuck{i}, rv]=qr(tuck{i}, 'lr');
    r1 = size(crx{i}, 1); r2 = size(crx{i}, 3); rt = size(crx{i},2);
    crx{i} = permute(crx{i}, [2, 1, 3]);
    crx{i} = reshape(crx{i}, rt, r1*r2);
    crx{i} = rv*crx{i};
    rt = size(rv,1);
    crx{i} = reshape(crx{i}, rt, r1, r2);
    crx{i} = permute(crx{i}, [2, 1, 3]);    
end;
% [crx,rv]=qr(crx, 'rl');
% crx = rv.'*crx;
for i=D:-1:2
    r1 = size(crx{i}, 1); rt = size(crx{i}, 2); r2 = size(crx{i}, 3);
    crx{i} = reshape(crx{i}, r1, rt*r2);
    [crx{i}, rv]=qr(crx{i}.', 0);
    r0 = size(crx{i-1}, 1); rt0 = size(crx{i-1}, 2);
    crx{i-1} = reshape(crx{i-1}, r0*rt0, r1);
    crx{i-1}=crx{i-1}*rv.';
    r1 = size(rv,1);
    crx{i-1} = reshape(crx{i-1}, r0, rt0, r1);
    crx{i} = reshape(crx{i}.', r1, rt, r2);
end;
for i=1:D
    r1 = size(crx{i}, 1); r2 = size(crx{i}, 3); rt = size(crx{i},2);
    crx{i} = permute(crx{i}, [1, 3, 2]);
    crx{i} = reshape(crx{i}, r1*r2, rt);
    [crx{i}, rv]=qr(crx{i}, 0);
    tuck{i} = tuck{i}*rv.';
    tuck{i} = qtt_fft1(tuck{i}, tol);
    % rank is reversed
    r = tuck{i}.r;
    rq1 = r(1);
    for k=1:rq1
        ek = [zeros(1,k-1), 1, zeros(1,rq1-k)];
        yk = ek*tuck{i}; % first rank is 1
        if (k==1)
            ynew = yk;
        else
            ynew = [ynew, yk];
            ynew = round(ynew, tol);
        end;
    end;
    [tuck{i}, rv]=qr(ynew, 'lr');
%     % reversing 2
%     % nonorth block is the last, but we need to SVD from the first
%     tuck{i} = core2cell(tuck{i});
%     for j=numel(tuck{i}):-1:2
%         tuck{i}{j} = reshape(tuck{i}{j}, r(j), 2*r(j+1));
%         [tuck{i}{j}, rv]=qr(tuck{i}{j}.', 0);
%         tuck{i}{j-1} = reshape(tuck{i}{j-1}, r(j-1)*2, r(j));
%         tuck{i}{j-1} = tuck{i}{j-1}*rv.';
%         r(j) = size(rv, 1);
%         tuck{i}{j} = reshape(tuck{i}{j}.', r(j), 2, r(j+1));
%         tuck{i}{j-1} = reshape(tuck{i}{j-1}, r(j-1), 2, r(j));
%     end;
%     tuck{i}{1} = permute(tuck{i}{1}, [2, 1, 3]);
%     tuck{i}{1} = reshape(tuck{i}{1}, 2, rq1*r(2));
%     [u,s,v]=svd(tuck{i}{1}, 'econ');
%     s = diag(s);
%     rnew = my_chop2(s, norm(s)*tol/sqrt(numel(tuck{i})-1));
%     u=u(:,1:rnew);
%     v=v(:,1:rnew)*diag(s(1:rnew));
%     tuck{i}{1} = reshape(u, 1, 2, rnew);
%     r(1) = 1;
%     rv = reshape(v', rnew*rq1, r(2));
%     for j=2:numel(tuck{i})
%         tuck{i}{j} = reshape(tuck{i}{j}, r(j), 2*r(j+1));
%         tuck{i}{j} = rv*tuck{i}{j};
%         r(j) = rnew;
%         tuck{i}{j} = reshape(tuck{i}{j}, r(j), rq1, 2, r(j+1));
%         tuck{i}{j} = permute(tuck{i}{j}, [1, 3, 2, 4]);
%         tuck{i}{j} = reshape(tuck{i}{j}, r(j)*2, rq1*r(j+1));
%         [u,s,v]=svd(tuck{i}{j}, 'econ');
%         s = diag(s);
%         rnew = my_chop2(s, norm(s)*tol/sqrt(numel(tuck{i})-1));
%         u=u(:,1:rnew);
%         v=v(:,1:rnew)*diag(s(1:rnew));
%         tuck{i}{j} = reshape(u, r(j), 2, rnew);
%         rv = reshape(v', rnew*rq1, r(j+1));        
%     end;
%     rv = reshape(rv, rnew, rq1);
%     tuck{i} = cell2core(tt_tensor, tuck{i});
    
    crx{i} = crx{i}*rv.';
    rt = size(rv,1);
    crx{i} = reshape(crx{i}, r1, r2, rt);
    crx{i} = permute(crx{i}, [1, 3, 2]);        
    if (i<D)
        crx{i} = reshape(crx{i}, r1*rt, r2);
        [crx{i}, rv]=qr(crx{i}, 0);
        r3 = size(crx{i+1}, 3); rt2 = size(crx{i+1},2);
        crx{i+1} = reshape(crx{i+1}, r2, rt2*r3);
        crx{i+1} = rv*crx{i+1};
        r2 = size(rv,1);
        crx{i} = reshape(crx{i}, r1, rt, r2);
        crx{i+1} = reshape(crx{i+1}, r2, rt2, r3);
    end;
end;

y = x;
y.core = cell2core(tt_tensor, crx);
y.tuck = tuck;
end