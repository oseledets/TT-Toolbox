function [y]=qtt_fft(x, ds, tol)
% Discrete Fourier Transform in the QTT format
%   function [y]=qtt_fft(x, ds, tol)
%   Computes the multidimensional DFT of the vector x in the QTT format with
%   the approximation tolerance tol.
%   DS denotes the splitting of the global tensor to the physical
%   variables, and should be the vector of quantics dimensions for each
%   physical mode.
%   Returns both the quantized and physical dimensions in the correct
%   order.
%
%   !!WARNING!! employs a (probably) costly rank reverse procedure.
%   
%   See S. Dolgov, B. Khoromskij, D. Savostyanov, 
%   Superfast Fourier transform using QTT approximation,
%   J. Fourier Anal. Appl., 18(5), 2012.


D = numel(ds);
y=[];
[x,nrm] = qr(x, 'rl');
x = nrm.'*x;
for i=1:D
    cury = chunk(x, 1, ds(i));
    if (i<D); x = chunk(x, ds(i)+1, x.d); end;
    
    cury = qtt_fft1(cury, tol);
    % tail ranks should be reversed
    r = cury.r;    
    r1 = r(1);
    r2 = r(ds(i)+1);    
    for k=1:r1
        ek = [zeros(1,k-1), 1, zeros(1,r1-k)];
        yk = ek*cury; % first rank is 1
        for m=1:r2
            em = [zeros(m-1,1); 1; zeros(r2-m,1)];
            if (m==1)
                ym = yk*em;
            else
                ym = [ym;(yk*em)]; % increasing the first rank
                ym = round(ym, tol);
            end;
        end;
        if (k==1)
            cury2 = ym;
        else
            cury2 = [cury2, ym]; % increasing the last rank
            cury2 = round(cury2, tol);
        end;
    end;
    [cury, rv]=qr(cury2, 'lr');    
%     % reversing 2
%     % nonorth block is the last
%     cury = core2cell(cury);    
%     cury{ds(i)} = permute(cury{ds(i)}, [1,3,2]);
%     cury{ds(i)} = reshape(cury{ds(i)}, r(ds(i))*r2, 2);
%     r(ds(i)+1)=1;
%     [u,s,v]=svd(cury{ds(i)}, 'econ');
%     s = diag(s);
%     rnew = my_chop2(s, norm(s)*tol/sqrt(ds(i)-1));
%     rv=u(:,1:rnew)*diag(s(1:rnew));
%     v=v(:,1:rnew);
%     cury{ds(i)} = reshape(v', rnew, 2);
%     rv = reshape(rv, r(ds(i)), r2*rnew);
%     for j=ds(i)-1:-1:2        
%         cury{j} = reshape(cury{j}, r(j)*2, r(j+1));
%         cury{j} = cury{j}*rv;
%         cury{j} = reshape(cury{j}, r(j), 2, r2, rnew);
%         cury{j} = permute(cury{j}, [1,3,2,4]);
%         r(j+1)=rnew;
%         cury{j} = reshape(cury{j}, r(j)*r2, 2*r(j+1));
%         [u,s,v]=svd(cury{j}, 'econ');
%         s = diag(s);
%         rnew = my_chop2(s, norm(s)*tol/sqrt(ds(i)-1));
%         rv=u(:,1:rnew)*diag(s(1:rnew));
%         v=v(:,1:rnew);
%         cury{j} = reshape(v', rnew, 2, r(j+1));
%         rv = reshape(rv, r(j), r2*rnew);
%     end;
%     cury{1} = reshape(cury{1}, r1*2, r(2));
%     cury{1} = cury{1}*rv;
%     r(2) = rnew;
%     cury{1} = reshape(cury{1}, r1, 2, r2, r(2));
%     cury{1} = permute(cury{1}, [3, 2, 1, 4]);
%     cury{1} = reshape(cury{1}, r2*2, r1*r(2));
%     [u,s,v]=svd(cury{1}, 'econ');
%     s = diag(s);
%     rnew = my_chop2(s, norm(s)*tol/sqrt(ds(i)-1));
%     u=u(:,1:rnew);
%     v=v(:,1:rnew)*diag(s(1:rnew));
%     cury{1} = reshape(u, r2, 2, rnew);
%     r(1) = r2;
%     rv = reshape(v', rnew*r1, r(2));
%     for j=2:ds(i)
%         cury{j} = reshape(cury{j}, r(j), 2*r(j+1));
%         cury{j} = rv*cury{j};
%         r(j) = rnew;
%         cury{j} = reshape(cury{j}, r(j), r1, 2, r(j+1));
%         cury{j} = permute(cury{j}, [1, 3, 2, 4]);
%         cury{j} = reshape(cury{j}, r(j)*2, r1*r(j+1));
%         [u,s,v]=svd(cury{j}, 'econ');
%         s = diag(s);
%         rnew = my_chop2(s, norm(s)*tol/sqrt(ds(i)-1));
%         u=u(:,1:rnew);
%         v=v(:,1:rnew)*diag(s(1:rnew));
%         cury{j} = reshape(u, r(j), 2, rnew);
%         rv = reshape(v', rnew*r1, r(j+1));        
%     end;
%     rv = reshape(rv, rnew, r1);
%     cury = cell2core(tt_tensor, cury);

    if (i<D)
        x=rv*x;
    else
        cury=cury*rv;
    end;
    y = kron(y, cury);
end;

end
