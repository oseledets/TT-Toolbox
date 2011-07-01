function [x] = dmrg_nlhad(B, M, A, f, x0, eps, rmax, nswp)

kickrank=2;

input_is_tt_tensor = 0;
if ( isa(f,'tt_tensor') )
  f=core(f);
  input_is_tt_tensor = 1;
  if ((nargin<5)||(isempty(x0)))
      x0=tt_random(tt_size(f), max(size(f)), 2);
  end;
else
    if ((nargin<5)||(isempty(x0)))
        x0=tt_random(tt_size(f), max(size(f)), 2);
    end;
end

if (isa(B,'tt_matrix') )
  B=core(B);
  input_is_tt_tensor = 1;
end
if ( isa(x0,'tt_tensor') )
  x0=core(x0);
  input_is_tt_tensor = 1;
end
if ( isa(A,'tt_matrix') )
  A=core(A);
  input_is_tt_tensor = 1;
end
if ( isa(M,'tt_matrix') )
  M=core(M);
  input_is_tt_tensor = 1;
end

if (nargin<8)||(isempty(nswp))
    nswp=1;
end;

x=x0;

x{1}=reshape(x{1}, size(x{1},1), 1, size(x{1},2));
f{1}=reshape(f{1}, size(f{1},1), 1, size(f{1},2));
A{1}=reshape(A{1}, size(A{1},1),size(A{1},2), 1, size(A{1},3)); %Bydlocode (@)
B{1}=reshape(B{1}, size(B{1},1),size(B{1},2), 1, size(B{1},3));
M{1}=reshape(M{1}, size(M{1},1),size(M{1},2), 1, size(M{1},3));

phB = cell(d+1,1); phB{1}=1; phB{d+1}=1;
phAM = cell(d+1,1); phAM{1}=1; phAM{d+1}=1;
phf = cell(d+1,1); phf{1}=1; phf{d+1}=1;

for swp=1:nswp
    % 1-to-d orthogonalization
    rvx = 1; rnewx=1; phBold=1; phAMold=1; phfold=1;

    for i=1:d-1
        cre = x{i};
        n1 = size(cre,1); rx1 = size(cre,2); rx2 = size(cre,3);
        cre = reshape(permute(cre, [2 1 3]), rx1, n1*rx2);
        cre = rvx*cre; % size rnew,n1,rx2
        rx1=rnewx;
        cre = reshape(cre, rx1*n1, rx2);
        [q,rvx]=qr(cre,0); % size rx1*n1,r2new - r2new,rx2
        rnewx = min(rx1*n1, rx2);
        x{i}=permute(reshape(q, rx1, n1, rnewx), [2 1 3]);
        
        % Phi for matrix B: phB{i}=phB{i-1} X^T B X
        b1 = B{i}; n1 = size(b1, 1); m1 = size(b1,2); rb1 = size(b1,3); rb2 = size(b1,4);
        x1 = x{i}; rx1 = size(x1, 2); rx2 = size(x1, 3);
        phBold = reshape(phBold, rx1*rb1, rx1);
        x1 = reshape(permute(x1, [2 1 3]), rx1, m1*rx2);
        phBold = phBold*x1; % size rx1'*rb1, m1*rx2
        phBold = reshape(phBold, rx1, rb1, m1, rx2);
        phBold = reshape(permute(phBold, [1 4 3 2]), rx1*rx2, m1*rb1);
        b1 = reshape(permute(b1, [2 3 1 4]), m1*rb1, n1*rb2);
        phBold = phBold*b1; % size rx1'*rx2, n1*rb2
        phBold = reshape(phBold, rx1, rx2, n1, rb2);
        phBold = reshape(permute(phBold, [1 3 4 2]), rx1*n1, rb2*rx2);
        x1 = reshape(x1, rx1*n1, rx2);
        phBold = (x1')*phBold; % size rx2', rb2*rx2
        phB{i+1} = reshape(phBold, rx2, rb2, rx2);
        
        % Phi for rhs: phf{i}=phf{i-1} X^T f
        f1 = f{i}; rf1 = size(f1, 2); rf2 = size(f1, 3);
        f1 = reshape(permute(f1, [2 1 3]), rf1, n1*rf2);
        phfold = phfold*f1; % size rx1', n1*rf2
        phfold = reshape(phfold, rx1*n1, rf2);
        phf{i+1} = (x1')*phfold; % size rx2, rf2
        
        % Now the hardest: very big phi for Mx .* Ax
        % phAMold = [rx1, rm1, ra1, rx1', rx1'']
        m1 = M{i}; rm1 = size(m1, 3); rm2 = size(m1, 4);
        a1 = A{i}; ra1 = size(a1, 3); ra2 = size(a1, 4);
        phAMold = reshape(phAMold, rx1*rm1*ra1*rx1, rx1);
        x1 = reshape(x1, rx1, m1*rx2);
        phAMold = phAMold*x1; % size rx1*rm1*ra1*rx1', m1*rx2''
        phAMold = reshape(phAMold, rx1, rm1, ra1, rx1, m1, rx2);
        phAMold = permute(phAMold, [1 2 4 6 5 3]);
        phAMold = reshape(phAMold, rx1*rm1*rx1*rx2, m1*ra1); % rx1*rm1*rx1'*rx2'',m1*ra1
        a1 = reshape(permute(a1, [2 3 1 4]), m1*ra1, n1*ra2);
        phAMold = phAMold*a1; % size rx1*rm1*rx1'*rx2'', n1*ra2
        phAMold = reshape(phAMold, rx1, rm1, rx1, rx2, n1, ra2);
        phAMold = permute(phAMold, [1 2 4 5 6 3]);
        phAMold = reshape(phAMold, rx1*rm1*rx2*n1*ra2, rx1);
        x1 = reshape(x1, rx1, m1*rx2);
        phAMold = phAMold*x1; % size rx1*rm1*rx2''*n1*ra2, m1'*rx2'
        phAMold = reshape(phAMold, rx1, rm1, rx2, n1, ra2, m1, rx2);
        phAMold = permute(phAMold, [1 3 4 5 7 6 2]);
        phAMold = reshape(phAMold, rx1*rx2*n1*ra2*rx2, m1*rm1); % rx1*rx2''*n1*ra2*rx2', m1'*rm1
        m1 = reshape(permute(m1, [2 3 1 4]), m1*rm1, n1*rm2);
        phAMold = phAMold*m1; % size rx1*rx2''*n1*ra2*rx2', n1'*rm2
        phAMold = reshape(phAMold, rx1, rx2, n1, ra2, rx2, n1, rm2);
        phAMold = permute(phAMold, [1 2 4 5 7 3 6]);
        phAMold = reshape(phAMold, rx1*rx2*ra2*rx2*rm2, n1, n1);
        % extract Hadamard products
        for j=1:n1
            phAMold(:,j,1) = phAMold(:,j,j);
        end;
        phAMold = reshape(phAMold(:,:,1), rx1, rx2, ra2, rx2, rm2, n1); % rx1*rx2''*ra2*rx2'*rm2*n1
        phAMold = permute(phAMold, [1 6 2 3 4 5]);
        phAMold = reshape(phAMold, rx1*n1, rx2*ra2*rx2*rm2);
        x1 = reshape(x1, rx1*n1, rx2);
        phAMold = (x1')*phAMold; % size rx2, rx2''*ra2*rx2'*rm2
        phAMold = reshape(phAMold, rx2, rx2, ra2, rx2, rm2);
        phAM{i+1} = permute(phAMold, [1 5 3 4 2]);
    end;
    
    
    % Now, start d-to-1 DMRG sweep
    phAMold=1; phfold=1; phBold = 1;
    for i=d:-1:2
        x1 = x{i-1}; m1 = size(x1,1); rx1 = size(x1,2); rx2 = size(x1,3);
        x2 = x{i}; m2 = size(x2,1); rx3 = size(x2,3);
        
        b1 = B{i-1}; n1 = size(b1,1); rb1 = size(b1,3); rb2 = size(b1,4);
        b2 = B{i}; n2 = size(b2,1); rb3 = size(b2,4);        
        % Projected B:
        locB = cell(2,1);
        locB{1} = phB{i-1};
        locB{1} = reshape(permute(locB{1}, [1 3 2]), rx1*rx1, rb1);
        b1 = reshape(permute(b1, [3 1 2 4]), rb1, n1*m1*rb2);
        locB{1} = locB{1}*b1; % size rx1'*rx1, n1*m1*rb2
        locB{1} = reshape(locB{1}, rx1, rx1, n1, m1, rb2);
        locB{1} = reshape(permute(locB{1}, [1 3 2 4 5]), rx1*n1, rx1*m1, rb2);
        
        locB{2} = phBold;
        locB{2} = reshape(permute(locB{2}, [1 3 2]), rx3*rx3, rb3);
        b2 = reshape(b2, n2*m2*rb2, rb3);
        locB{2} = locB{2}*(b2.'); % size rx3'*rx3, n2*m2*rb2
        locB{2} = reshape(locB{2}, rx3, rx3, n2, m2, rb2);
        locB{2} = reshape(permute(locB{2}, [3 1 4 2 5]), n2*rx3, m2*rx3, rb2);        
        
        % Projected rhs
        rhs = phf{i-1};
        f1 = f{i-1}; rf1 = size(f1,2); rf2 = size(f1,3);
        f2 = f{i}; rf3 = size(f2,3);        
        
        f1 = reshape(permute(f1, [2 1 3]), rf1, n1*rf2);
        rhs = rhs*f1; % size rx1, n1*rf2
        rhs = reshape(rhs, rx1*n1, rf2);
        f2 = reshape(permute(f2, [2 1 3]), rf2, n2*rf3);
        rhs = rhs*f2; % size rx1*n1, n2*rf3
        rhs = reshape(rhs, rx1*n1*n2, rf3);
        rhs = rhs*(phfold.'); % size rx1*n1*n2, rx3
        rhs = reshape(rhs, rx1*n1*n2*rx3, 1);
        
        % projected nonlinear 2-3D tensor
        a1 = A{i-1};  ra1 = size(a1,3); ra2 = size(a1,4);
        a2 = A{i};  ra3 = size(a2,4);                        
        m1 = M{i-1};  rm1 = size(m1,3); rm2 = size(m1,4);
        m2 = M{i};  rm3 = size(m2,4);                        
        
        AM = cell(2,1);
        AM{1} = phAM{i-1}; % rx1, rm1, ra1, rx1'(m), rx1''(a)
        AM{1} = reshape(permute(AM{1}, [1 4 5 2 3]), rx1*rx1*rx1*rm1, ra1);
        a1 = reshape(permute(a1, [3 1 2 4]), ra1, n1*m1*ra2);
        AM{1}=AM{1}*a1; % size rx1*rx1'*rx1''*rm1, n1*m1*ra2
        AM{1}=reshape(AM{1}, rx1, rx1, rx1, rm1, n1, m1, ra2);
        AM{1}=reshape(permute(AM{1}, [1 2 3 5 6 7 4]), rx1*rx1*rx1*n1*m1*ra2, rm1);
        m1 = reshape(permute(m1, [3 1 2 4]), rm1, n1*m1*rm2);
        AM{1}=AM{1}*m1; % size rx1*rx1'*rx1''*n1*m1*ra2, n1'*m1'*rm2
        AM{1}=reshape(AM{1}, rx1*rx1*rx1, n1, m1, ra2, n1, m1, rm2);
        AM{1}=permute(AM{1}, [1 4 7 3 6 2 5]);
        % rx1*rx1'*rx1''*ra2*rm2*m1(a)*m1'(m), n1, n1
        AM{1}=reshape(AM{1}, rx1*rx1*rx1*ra2*rm2*m1*m1, n1, n1);
        for j=1:n1
            AM{1}(:,j,1)=AM{1}(:,j,j);
        end;
        AM{1}=reshape(AM{1}(:,:,1), rx1,rx1,rx1,ra2,rm2,m1,m1, n1);
        AM{1}=permute(AM{1}, [[1 8] [2 7] [3 6] [4 5]]);
        AM{1}=reshape(AM{1}, rx1*n1, rx1*m1, rx1*m1, ra2*rm2);
        
        AM{2} = phAMold; % rx3, rm3, ra3, rx3'(m), rx3''(a)
        AM{2} = reshape(permute(AM{2}, [1 4 5 2 3]), rx3*rx3*rx3*rm3, ra3);
        a2 = reshape(a2, n2*m2*ra2, ra3);
        AM{2}=AM{2}*(a2.'); % size rx3*rx3'*rx3''*rm3, n2*m2*ra2
        AM{2}=reshape(AM{2}, rx3, rx3, rx3, rm3, n2, m2, ra2);
        AM{2}=reshape(permute(AM{2}, [1 2 3 5 6 7 4]), rx3*rx3*rx3*n2*m2*ra2, rm3);
        m2 = reshape(m2, n2*m2*rm2, rm3);
        AM{2}=AM{2}*(m2.'); % size rx3*rx3'*rx3''*n2*m2*ra2, n2'*m2'*rm2
        AM{2}=reshape(AM{2}, rx3*rx3*rx3, n2, m2, ra2, n2, m2, rm2);
        AM{2}=permute(AM{2}, [1 4 7 3 6 2 5]);
        % rx3*rx3'*rx3''*ra2*rm2*m2(a)*m2'(m), n2, n2
        AM{2}=reshape(AM{2}, rx3*rx3*rx3*ra2*rm2*m2*m2, n2, n2);
        for j=1:n2
            AM{2}(:,j,1)=AM{2}(:,j,j);
        end;
        AM{2}=reshape(AM{2}(:,:,1), rx3,rx3,rx3,ra2,rm2,m2,m2, n2);
        AM{2}=permute(AM{2}, [[8 1] [7 2] [6 3] [4 5]]);
        AM{2}=reshape(AM{2}, n2*rx3, m2*rx3, m2*rx3, ra2*rm2);  
        
        
        % Put newton minimizer here
        %
        %
        %
        sol = reshape(sol, rx1*m1, m2*rx3);
        [u,s,v]=svd(sol, 'econ');
        r = my_chop2(diag(s), norm(diag(s))*eps);
        r = min(r, rmax);
        v = [v(:,1:r), randn(m2*rx3, kickrank)];
        u = u(:,1:r)*s(1:r,1:r);
        u = [u zeros(rx1*m, kickrank)];
        [v,rv]=qr(v, 0);
        r = size(v,2);
        u=u*(rv.');
        x{i}=permute(reshape(v, m2, rx3, r), [1 3 2]);
        x{i-1}=permute(reshape(u, rx1, m1, r), [2 1 3]);
        
        
        % Update phi
        % Phi for matrix B: phB{i}=phB{i+1} X^T B X
        b1 = B{i}; n1 = size(b1, 1); m1 = size(b1,2); rb1 = size(b1,3); rb2 = size(b1,4);
        x1 = x{i}; rx1 = size(x1, 2); rx2 = size(x1, 3);
        phBold = reshape(phBold, rx2*rb2, rx2);
        x1 = reshape(permute(x1, [2 1 3]), rx1*m1, rx2);
        phBold = phBold*(x1.'); % size rx2'*rb2, m1*rx1
        phBold = reshape(phBold, rx2, rb2, m1, rx1);
        phBold = reshape(permute(phBold, [1 4 3 2]), rx2*rx1, m1*rb2);
        b1 = reshape(permute(b1, [2 4 1 3]), m1*rb2, n1*rb1);
        phBold = phBold*b1; % size rx2'*rx1, n1*rb1
        phBold = reshape(phBold, rx2, rx1, n1, rb1);
        phBold = reshape(permute(phBold, [3 1 4 2]), n1*rx2, rb1*rx1);
        x1 = reshape(x1, rx1, n1*rx2);
        phBold = conj(x1)*phBold; % size rx1', rb1*rx1
        phB{i} = reshape(phBold, rx1, rb1, rx1);
        
        % Phi for rhs: phf{i}=phf{i+1} X^T f
        f1 = f{i}; rf1 = size(f1, 2); rf2 = size(f1, 3);
        f1 = reshape(f1, n1*rf1, rf2);
        phfold = phfold*(f1.'); % size rx2', n1*rf1
        phfold = reshape(phfold, rx2, n1, rf1);
        phfold = permute(phfold, [2 1 3]);
        phfold = reshape(phfold, n1*rx2, rf1);
        phf{i} = conj(x1)*phfold; % size rx1, rf1
        
        % Now the hardest: very big phi for Mx .* Ax
        % phAMold = [rx2, rm2, ra2, rx2', rx2'']
        m1 = M{i}; rm1 = size(m1, 3); rm2 = size(m1, 4);
        a1 = A{i}; ra1 = size(a1, 3); ra2 = size(a1, 4);
        phAMold = reshape(phAMold, rx2*rm2*ra2*rx2, rx2);
        x1 = reshape(x1, rx1*m1, rx2);
        phAMold = phAMold*(x1.'); % size rx2*rm2*ra2*rx2', rx1''*m1
        phAMold = reshape(phAMold, rx2, rm2, ra2, rx2, rx1, m1);
        phAMold = permute(phAMold, [1 2 4 5 6 3]);
        phAMold = reshape(phAMold, rx2*rm2*rx2*rx1, m1*ra2); % rx2*rm2*rx2'*rx1'',m1*ra2
        a1 = reshape(permute(a1, [2 4 1 3]), m1*ra2, n1*ra1);
        phAMold = phAMold*a1; % size rx2*rm2*rx2'*rx1'', n1*ra1
        phAMold = reshape(phAMold, rx2, rm2, rx2, rx1, n1, ra1);
        phAMold = permute(phAMold, [1 2 4 5 6 3]);
        phAMold = reshape(phAMold, rx2*rm2*rx1*n1*ra1, rx2);
        x1 = reshape(x1, rx1*m1, rx2);
        phAMold = phAMold*(x1.'); % size rx2*rm2*rx1''*n1*ra1, rx1'*m1'
        phAMold = reshape(phAMold, rx2, rm2, rx1, n1, ra1, rx1, m1);
        phAMold = permute(phAMold, [1 3 4 5 6 7 2]);
        phAMold = reshape(phAMold, rx2*rx1*n1*ra1*rx1, m1*rm2); % rx2*rx1''*n1*ra1*rx1', m1'*rm2
        m1 = reshape(permute(m1, [2 4 1 3]), m1*rm2, n1*rm1);
        phAMold = phAMold*m1; % size rx2*rx1''*n1*ra1*rx1', n1'*rm1
        phAMold = reshape(phAMold, rx2, rx1, n1, ra1, rx1, n1, rm1);
        phAMold = permute(phAMold, [1 2 4 5 7 3 6]);
        phAMold = reshape(phAMold, rx2*rx1*ra1*rx1*rm1, n1, n1);
        % extract Hadamard products
        for j=1:n1
            phAMold(:,j,1) = phAMold(:,j,j);
        end;
        phAMold = reshape(phAMold(:,:,1), rx2, rx1, ra1, rx1, rm1, n1); % rx2*rx1''*ra1*rx1'*rm1*n1
        phAMold = permute(phAMold, [6 1 2 3 4 5]);
        phAMold = reshape(phAMold, n1*rx2, rx1*ra1*rx1*rm1);
        x1 = reshape(x1, rx1, n1*rx2);
        phAMold = conj(x1)*phAMold; % size rx1, rx1''*ra1*rx1'*rm1
        phAMold = reshape(phAMold, rx1, rx1, ra1, rx1, rm1);
        phAM{i} = permute(phAMold, [1 5 3 4 2]);        
    end;
end;

end