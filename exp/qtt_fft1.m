function [y]=qtt_fft1(x, tol)
% Discrete Fourier Transform in the QTT format, 1D
%   function [y]=qtt_fft1(x, tol)
%   Computes the one-dimensional DFT of the vector x in the QTT format with
%   the approximation tolerance tol.
%   Returns the frequency bits in the correct ordering. The tail ranks
%   (if ~=1) will be reversed.
%   
%   See Dolgov,Khoromskij,Savostyanov, 
%   "Superfast Fourier Transform Using QTT Approximation" (2012)
%   for further details.

d = x.d;
r = x.r;
y = core2cell(x);
for i=d:-1:2
    r1 = size(y{i},1); r2 = size(y{i},3);
    crd2 = zeros(r1, 2, r2);
    % last block +-
    crd2(:,1,:)=(y{i}(:,1,:)+y{i}(:,2,:))/sqrt(2);
    crd2(:,2,:)=(y{i}(:,1,:)-y{i}(:,2,:))/sqrt(2);    
    % last block twiddles
    y{i} = zeros(r1*2,2,r2);
    y{i}(1:r1, 1, 1:r2) = crd2(:,1,:);
    y{i}(r1+1:r1*2, 2, 1:r2) = crd2(:,2,:);
    % 1..i-1 block twiddles and qr
    rv = 1;
    for j=1:i-1
        cr = y{j};        
        r1 = size(cr,1); r2 = size(cr,3);
        if (j==1)
            r(j)=r1; r(j+1)=r2*2;
            y{j} = zeros(r(j), 2, r(j+1));
            y{j}(1:r1, :, 1:r2) = cr;
            y{j}(1:r1, 1, r2+1:r(j+1)) = cr(:,1,:);
            y{j}(1:r1, 2, r2+1:r(j+1)) = exp(-2*pi*1i/(2^(i-j+1)))*cr(:,2,:);
        else    
            r(j)=r1*2; r(j+1)=r2*2;
            y{j} = zeros(r(j), 2, r(j+1));
            y{j}(1:r1, :, 1:r2) = cr;
            y{j}(r1+1:r(j), 1, r2+1:r(j+1)) = cr(:,1,:);
            y{j}(r1+1:r(j), 2, r2+1:r(j+1)) = exp(-2*pi*1i/(2^(i-j+1)))*cr(:,2,:);            
        end;        
        y{j} = reshape(y{j}, r(j), 2*r(j+1));
        y{j} = rv*y{j};
        r(j) = size(y{j}, 1);
        y{j} = reshape(y{j}, r(j)*2, r(j+1));
        [y{j}, rv]=qr(y{j}, 0);      
        y{j} = reshape(y{j}, r(j), 2, size(rv, 1));
    end;
    y{i} = reshape(y{i}, r(i), 2*r(i+1));
    y{i} = rv*y{i};
    r(i) = size(rv,1);
    % backward svd
    for j=i:-1:2
        [u,s,v]=svd(y{j}, 'econ');
        s = diag(s);
        rnew = my_chop2(s, norm(s)*tol/sqrt(i-1));
        u = u(:,1:rnew)*diag(s(1:rnew));
        v=v(:,1:rnew);
        y{j} = reshape(v', rnew, 2, r(j+1));
        y{j-1} = reshape(y{j-1}, r(j-1)*2, r(j));
        y{j-1}=y{j-1}*u;
        r(j) = rnew;
        y{j-1} = reshape(y{j-1}, r(j-1), 2*r(j));
    end;    
    y{1} = reshape(y{1}, r(1), 2, r(2));
end;
% FFT on the first block
y{1} = permute(y{1}, [2,1,3]);
y{1} = reshape(y{1}, 2, r(1)*r(2));
y{1} = [1,1;1,-1]*y{1}/sqrt(2);
y{1} = reshape(y{1}, 2, r(1), r(2));
y{1} = permute(y{1}, [2,1,3]);

% Reverse the train
y2 = cell(d,1);
for i=1:d
    y2{d-i+1} = permute(y{i}, [3,2,1]);
end;

y = cell2core(tt_tensor, y2);
end