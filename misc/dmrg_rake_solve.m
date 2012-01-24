function [xf, xc]=dmrg_rake_solve(Af, Ac, yf, yc, tol, varargin)

nswp = 10;
local_format = 'full';
% local_format = 'tt';
max_full_size = 1500;
nrestart = 40;
gmres_iters = 2;
verb = 2;
kickrank = 1;

d = yc.d; % Physical dim.
L = zeros(1,d); % Quantics dims
n = zeros(max(L), d); % Physical mode sizes
for i=1:d
    L(i) = yf{i}.d;
    n(1:L(i), i) = yf{i}.n;
end;

xc = tt_rand(10,d,10);
xf = cell(d,1);
for i=1:d
    xf{i} = tt_rand(n(1:L(i),i), L(i), [1;10*ones(L(i),1)]);
end;

% Extract ranks; Note that rf(L(i)+1,i) = r_tuck(i)
rcy = yc.r;
rfy = zeros(max(L)+1, d);
for i=1:d
    rfy(1:L(i)+1, i) = yf{i}.r;
end;
rcA = Ac.r;
rfA = zeros(max(L)+1, d);
for i=1:d
    rfA(1:L(i)+1, i) = Af{i}.r;
end;
rcx = xc.r;
rfx = zeros(max(L)+1, d);
for i=1:d
    rfx(1:L(i)+1, i) = xf{i}.r;
end;

% Init phis. Thousands of them... (c)
phcAl = cell(d+1,1); phcAl{1}=1; % Core, matrix, left
phcAr = cell(d+1,1); phcAr{d+1}=1; % Core, matrix, right
phfAb = cell(d+1,1); % Factors, matrix, bottom
for i=1:d
    phfAb{i} = cell(L(i)+1,1);
    phfAb{i}{1}=1;
end;
phfAt = cell(d+1,1); % Factors, matrix, top
for i=1:d
    phfAt{i} = cell(L(i)+1,1);
    % Very bad: L+1-th block of "top" phi has tucker ranks instead of 1
end;
phcyl = cell(d+1,1); phcyl{1}=1; % Core, rhs, left
phcyr = cell(d+1,1); phcyr{d+1}=1; % Core, rhs, right
phfyb = cell(d+1,1); % Factors, rhs, bottom
for i=1:d
    phfyb{i} = cell(L(i)+1,1);
    phfyb{i}{1}=1;
end;
phfyt = cell(d+1,1); % Factors, rhs, top
for i=1:d
    phfyt{i} = cell(L(i)+1,1);
    % Very bad: L+1-th block of "top" phi has tucker ranks instead of 1
end;

for swp=1:nswp
    dx_max = 0;
    res_max = 0;
    r_max = 0;
%     xc = xc + tt_rand(xc.n, d, 2);
%     rcx = xc.r;
%     for i=1:d
%         curcr = xf{i}{L(i)};
%         curcr(:,:,rfx(L(i)+1,i)+1:rfx(L(i)+1,i)+2) = randn(rfx(L(i),i), n(L(i),i), 2);
%         xf{i}{L(i)} = curcr;
%         curcr = xc{i};
%         curcr(:,rfx(L(i)+1,i)+1:rfx(L(i)+1,i)+2,:) = randn(rcx(i), 2, rcx(i+1));
%         xc{i} = curcr;
%         rfx(1:L(i)+1,i)=xf{i}.r;
%     end;
%     rcx = xc.r;
    % Right-to-left, bottom-to-top QR and phis
    for i=d:-1:1 % physical dims/core
        for j=1:L(i) % quantics dims
            cr = xf{i}{j};
            cr = reshape(cr, rfx(j,i)*n(j,i), rfx(j+1,i));
            [cr, rv] = qr(cr, 0);
            % What is our next core?
            if (j<L(i))
                % we are still on a "tooth"
                cr2 = xf{i}{j+1};
                cr2 = reshape(cr2, rfx(j+1,i), n(j+1,i)*rfx(j+2,i));
                cr2 = rv*cr2;
                rfx(j+1,i) = size(cr, 2);
                xf{i}{j} = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
                xf{i}{j+1} = reshape(cr2, rfx(j+1,i), n(j+1,i), rfx(j+2,i));
            else
                % We have to convlove rv to the tucker core
                cr2 = xc{i};
                cr2 = permute(cr2, [2, 1, 3]);
                cr2 = reshape(cr2, rfx(j+1,i), rcx(i)*rcx(i+1));
                cr2 = rv*cr2;
                rfx(j+1,i) = size(cr, 2);
                cr2 = reshape(cr2, rfx(j+1,i), rcx(i), rcx(i+1));
                cr2 = permute(cr2, [2, 1, 3]);
                xf{i}{j} = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
                xc{i} = cr2;
            end;
            % Update bottom phis
            cr = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
            phfAb{i}{j+1} = compute_next_Phi(phfAb{i}{j}, cr, Af{i}{j}, cr, 'lr');
%             curph = phfAb{i}{j};
%             curph = reshape(curph, rfx(j,i)*rfA(j,i), rfx(j,i));
%             cr = reshape(cr, rfx(j,i), n(j,i)*rfx(j+1,i));
%             curph = curph*cr; % size rfx1*rfA1, n*rfx2
%             a1 = Af{i}{j};
%             a1 = permute(a1, [2, 4, 1, 3]);
%             a1 = reshape(a1, n(j,i)*rfA(j+1,i), rfA(j,i)*n(j,i));
%             curph = reshape(curph, rfx(j,i), rfA(j,i)*n(j,i), rfx(j+1,i));
%             curph = permute(curph, [2, 1, 3]);
%             curph = reshape(curph, rfA(j,i)*n(j,i), rfx(j,i)*rfx(j+1,i));
%             curph = a1*curph; % size n*rfA2, rfx1*rfx2
%             curph = reshape(curph, n(j,i), rfA(j+1,i), rfx(j,i), rfx(j+1,i));
%             curph = permute(curph, [3, 1, 2, 4]);
%             curph = reshape(curph, rfx(j,i)*n(j,i), rfA(j+1,i)*rfx(j+1,i));
%             cr = reshape(cr, rfx(j,i)*n(j,i), rfx(j+1,i));
%             curph = (cr')*curph;
%             phfAb{i}{j+1} = reshape(curph, rfx(j+1,i), rfA(j+1,i), rfx(j+1,i));
            
            phfyb{i}{j+1} = compute_next_Phi(phfyb{i}{j}, cr, [], yf{i}{j}, 'lr');
%             curph = phfyb{i}{j};
%             y1 = yf{i}{j};            
%             y1 = reshape(y1, rfy(j,i), n(j,i)*rfy(j+1,i));
%             curph = curph*y1; % size rfx1, n*rfy2
%             curph = reshape(curph, rfx(j,i)*n(j,i), rfy(j+1,i));
%             curph = (cr')*curph;
%             phfyb{i}{j+1} = curph;
        end;
        if (i>1) % QR to the next tucker core
            % cr2 is what we need
            % n here is in fact a tucker rank
            cr2 = reshape(cr2, rcx(i), rfx(L(i)+1,i)*rcx(i+1));
            [cr2, rv]=qr(cr2.', 0);
            cr3 = xc{i-1};
            cr3 = reshape(cr3, rcx(i-1)*rfx(L(i-1)+1,i-1), rcx(i));
            cr3 = cr3*(rv.');
            rcx(i) = size(cr2, 2);
            cr2 = cr2.'; % sizes rcx1, rtuck*rcx2
            xc{i} = reshape(cr2, rcx(i), rfx(L(i)+1,i), rcx(i+1));
            xc{i-1} = reshape(cr3, rcx(i-1), rfx(L(i-1)+1,i-1), rcx(i));
            
            % Update right phis
            cr2 = reshape(cr2, rcx(i), rfx(L(i)+1,i), rcx(i+1));
            a1 = core_matrix(Ac{i}, phfAb{i}{L(i)+1}); % size rcA1, rtx, rtx, rcA2
            phcAr{i} = compute_next_Phi(phcAr{i+1}, cr2, a1, cr2, 'rl');
            % New size: rcx1, rcA1, rcx1
            
            y1 = core_vector(yc{i}, phfyb{i}{L(i)+1}); % size rcy1, rtx, rcy2
            phcyr{i} = compute_next_Phi(phcyr{i+1}, cr2, [], y1, 'rl');
            % New size: rcx1, rcy1
            
%             curph = phcAr{i+1};
%             curph = reshape(curph, rcx(i+1)*rcA(i+1), rcx(i+1));
%             cr2 = reshape(cr2, rcx(i)*rtx, rcx(i+1));
%             curph = curph*(cr.'); % size rcx2*rcA2, rcx1*rtuckx
%             curph = reshape(curph, rcx(i+1), rcA(i+1), rcx(i), rtx);
%             curph = permute(curph, [4, 2, 3, 1]);
%             curph = reshape(curph, rtx*rcA(i+1), rcx(i)*rcx(i+1));
%             a1 = Ac{i}; % sizes rca1, rtucka, rca2
%             a1 = permute(a1, [2, 1, 3]);
%             a1 = reshape(a1, rfA(L(i)+1,i), rcA(i)*rcA(i+1));
%             phbot = phfAb{i}{L(i)+1}; % sizes rtuckx, rtucka, rtuckx
%             phbot = permute(phbot, [1, 3, 2]);
%             phbot = reshape(phbot, rtx*rtx, rfA(L(i)+1,i));
%             a1 = phbot*a1; % size rtuckx'*rtuckx, rca1*rca2
%             a1 = reshape(a1, rtx, rtx, rcA(i), rcA(i+1));
%             a1 = permute(a1, [1, 3, 2, 4]);
%             a1 = reshape(a1, rtx*rcA(i), rtx*rcA(i+1));
%             curph = a1*curph; % size rtx*rcA1, rcx1*rcx2
%             curph = reshape(curph, rtx, rcA(i), rcx(i), rcx(i+1));
%             curph = permute(curph, [1, 4, 2, 3]);
%             curph = reshape(curph, rtx*rcx(i+1), rcA(i)*rcx(i));
%             cr2 = reshape(cr2, rcx(i), rtx*rcx(i+1));
%             curph = conj(cr2)*curph;
%             phcAr{i} = reshape(curph, rcx(i), rcA(i), rcx(i));
%            
%             
%             curph = phcyr{i+1};
%             y1 = yf{i}{j};            
%             y1 = reshape(y1, rfy(j,i), n(j,i)*rfy(j+1,i));
%             curph = curph*y1; % size rfx1, n*rfy2
%             curph = reshape(curph, rfx(j,i)*n(j,i), rfy(j+1,i));
%             curph = (cr')*curph;
%             phfAb{i}{j+1} = curph;            
        end;
    end;
    
    % Now we are at the position {1,L(1)+1}
    for i=1:d
        % DMRG over the factor
        % Convolve the core block to the last factor block
        crb = xc{i};
        crb = permute(crb, [2, 1, 3]);
        crb = reshape(crb, rfx(L(i)+1,i), rcx(i)*rcx(i+1));
        crf = xf{i}{L(i)};
        crf = reshape(crf, rfx(L(i),i)*n(L(i),i), rfx(L(i)+1,i));
        crf = crf*crb;
        crf = reshape(crf, rfx(L(i),i), n(L(i),i), rcx(i)*rcx(i+1));
        xf{i}{L(i)} = crf;
        
        for j=L(i):-1:2
            % Rename some quantities for brevity
            rx1 = rfx(j-1,i); ra1 = rfA(j-1,i); ry1 = rfy(j-1,i); n1 = n(j-1,i);            
            n2 = n(j,i);
            % if we are touching the L(i)+1-th core = tucker core
            if (j==L(i))
                rx2 = rfx(L(i),i);
                ra2 = rfA(L(i)+1,i); 
                rx3 = rcx(i)*rcx(i+1); 

                % Preparing the matrix and rhs in this case smells like
                % manure.
                % phcAl - Ac{i} - phcAr     <- this has to be convolved
                %          |                <- this has to be convolved
                %         Af{i}{j}
                %          |
                %         Af{i}{j-1}
                %          |                <- this has to be convolved
                %         phfAb{j-1}
                lefta1 = permute(phcAl{i}, [1, 3, 2]);
                lefta1 = reshape(lefta1, rcx(i)*rcx(i), rcA(i)); % rcx1'*rcx1, rcA1
                lefta1 = lefta1*reshape(Ac{i}, rcA(i), ra2*rcA(i+1));
                % We will need it later to compute a left local block for
                % tucker
                lefta1 = reshape(lefta1, rcx(i)*rcx(i)*ra2, rcA(i+1));
                ph2 = permute(phcAr{i+1}, [2, 1, 3]);
                ph2 = reshape(ph2, rcA(i+1), rx3*rx3);
                a1 = lefta1*ph2; % size rcx(i)*rcx(i)*rfA(L(i)+1,i), rcx(i+1)*rcx(i+1))
                a1 = reshape(a1, n2, n2, ra2, rx3, rx3);
                a1 = permute(a1, [1, 4, 2, 5, 3]);
                % First block of the local matrix is ready
                a1 = reshape(a1, n2*rx3, n2*rx3, ra2);
                
                lefty1 = phcyl{i};
                lefty1 = lefty1*reshape(yc{i}, rcy(i), ry2*rcy(i+1));
                % We will need it later to compute left tucker block
                lefty1 = reshape(lefty1, rcx(i)*ry2, rcy(i+1));
                ph2 = phcyr{i+1}.';
                y1 = lefty1*ph2; % size rcx1*rty, rcx2
                y1 = reshape(y1, rcx(i), ry2, rcx(i+1));
                y1 = permute(y1, [1, 3, 2]);
                % First block of the rhs is ready
                y1 = reshape(y1, n2*rx3, ry2);
                
                x1 = xc{i};
                x1 = permute(x1, [1, 3, 2]);
                x1 = reshape(x1, n2*rx3, rx2);
            else
                % We are inside the factor
                n2 = n(j,i); rx2 = rfx(j,i);
                ra2 = rfA(j,i); ry2 = rfy(j,i);
                rx3 = rfx(j+1,i); ra3 = rfA(j+1,i); ry3 = rfy(j+1,i);
                
                a1 = permute(phfAt{i}{j+1}, [2, 1, 3]);
                a1 = reshape(a1, ra3, rx3*rx3); 
                a1 = reshape(Af{i}{j}, ra2*n2*n2, ra3)*a1;
                % We will need it later to compute next phfAt
                a1 = reshape(a1, ra2, n2, n2, rx3, rx3);
                a1 = permute(a1, [2, 4, 3, 5, 1]);
                % First block of the local matrix is ready
                a1 = reshape(a1, n2*rx3, n2*rx3, ra2);
                
                y1 = phfyt{i}{j+1}.';
                y1 = reshape(yf{i}{j}, ry2*n2, ry3)*y1;
                % We will need it later to compute next phfyt
                y1 = reshape(y1, ry2, n2*rx3);
                % First block of the rhs is ready
                y1 = y1.';
                
                x1 = xf{i}{j};
                x1 = reshape(x1, rx2, n2*rx3);                
                x1 = x1.';
            end;
            
            a2 = phfAb{i}{j-1};
            a2 = permute(a2, [1, 3, 2]);
            a2 = reshape(a2, rx1*rx1, ra1);
            a2 = a2*reshape(Af{i}{j-1}, ra1, n1*n1*ra2); % rx1'*rx1,n1'*n1,rta
            a2 = reshape(a2, rx1, rx1, n1, n1, ra2);
            a2 = permute(a2, [1, 3, 2, 4, 5]);
            % Second block of the local matrix is ready
            a2 = reshape(a2, rx1*n1, rx1*n1, ra2);
            
            y2 = phfyb{i}{j-1};
            y2 = y2*reshape(yf{i}{j-1}, ry1, n1*ry2); % rx1'*rx1,n1'*n1,rta
            % Second block of the rhs is ready
            y2 = reshape(y2, rx1*n1, ry2);
            
            x2 = xf{i}{j-1};
            x2 = reshape(x2, rx1*n1, rx2);
            
            if (verb>1)
                fprintf('=factor= swp %d, block {%d}{%d}, ', swp, i, j);
            end;
            [u,s,v,r,dx_max,res_max]=local_solve(a2, a1, y2, y1, x2, x1, ...
                rx1, n1, n2, rx3, ra2, ...
                tol / sqrt(sum(L)) / 2, res_max, dx_max, ...
                local_format, max_full_size, nrestart, gmres_iters, verb);
            u = u*s;
            
            % kick
            v = reort(v, randn(n2*rx3, kickrank));
            radd = size(v,2)-r;
            u = [u, zeros(rx1*n1, radd)];
            r = size(v,2);
            
            r_max = max(r_max, r);

            % Recompute phis
            % Parts of phis are ready in a1, y1
            a1 = permute(a1, [1, 3, 2]);
            a1 = reshape(a1, n2*rx3*ra2, n2*rx3);
            curph = a1*v;
            curph = reshape(curph, n2*rx3, ra2*r);
            curph = (v')*curph;
            phfAt{i}{j} = reshape(curph, r, ra2, r);            
            phfyt{i}{j} = (v')*y1;
            
            % Stuff back
            rfx(j,i) = r;
            xf{i}{j-1} = reshape(u, rx1, n1, r);
            v = reshape(v, n2, rx3, r);
            if (j==L(i)+1)
                % v is the tucker core, u is the last factor's
                xc{i} = permute(v, [1, 3, 2]);                
            else
                xf{i}{j} = permute(v, [3, 1, 2]);
            end;
        end;
        
        % Go back to the core, QR @ Phi!
        for j=1:L(i) % quantics dims
            cr = xf{i}{j};
            cr = reshape(cr, rfx(j,i)*n(j,i), rfx(j+1,i));
            [cr, rv] = qr(cr, 0);
            % What is our next core?
            if (j<L(i))
                % we are still on a "tooth"
                cr2 = xf{i}{j+1};
                cr2 = reshape(cr2, rfx(j+1,i), n(j+1,i)*rfx(j+2,i));
                cr2 = rv*cr2;
                rfx(j+1,i) = size(cr, 2);
                xf{i}{j} = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
                xf{i}{j+1} = reshape(cr2, rfx(j+1,i), n(j+1,i), rfx(j+2,i));
            else
                % We have to convlove rv to the tucker core
                cr2 = xc{i};
                cr2 = permute(cr2, [2, 1, 3]);
                cr2 = reshape(cr2, rfx(j+1,i), rcx(i)*rcx(i+1));
                cr2 = rv*cr2;
                rfx(j+1,i) = size(cr, 2);
                cr2 = reshape(cr2, rfx(j+1,i), rcx(i), rcx(i+1));
                cr2 = permute(cr2, [2, 1, 3]);
                xf{i}{j} = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
                xc{i} = cr2;
            end;
            % Update bottom phis
            cr = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
            phfAb{i}{j+1} = compute_next_Phi(phfAb{i}{j}, cr, Af{i}{j}, cr, 'lr');
            phfyb{i}{j+1} = compute_next_Phi(phfyb{i}{j}, cr, [], yf{i}{j}, 'lr');
        end;
        
        % DMRG over tucker's edge
        if (i<d)
            rx1 = rcx(i); rx2 = rcx(i+1); rx3 = rcx(i+2);
            n1 = rfx(L(i)+1, i); n2 = rfx(L(i+1)+1,i+1);
            rta1 = rfA(L(i)+1, i); rta2 = rfA(L(i+1)+1,i+1);
            rty1 = rfy(L(i)+1, i); rty2 = rfy(L(i+1)+1,i+1);
            
            % We have lefta1 from the left, and phfAb from the bottom
            a1 = reshape(lefta1, rx1*rx1, rta1, rcA(i+1));
            a1 = permute(a1, [1, 3, 2]);
            a1 = reshape(a1, rx1*rx1*rcA(i+1), rta1);
            ph2 = phfAb{i}{L(i)+1};
            ph2 = permute(ph2, [2, 1, 3]);
            ph2 = reshape(ph2, rta1, n1*n1);
            a1 = a1*ph2;
            a1 = reshape(a1, rx1, rx1, rcA(i+1), n1, n1);
            a1 = permute(a1, [1, 4, 2, 5, 3]);
            % First matrix block is ready
            a1 = reshape(a1, rx1*n1, rx1*n1, rcA(i+1));
            
            % For the second block:  -A2 - phcAr=
            %                          |
            %                        -phfAb-
            a2 = phcAr{i+2};
            a2 = permute(a2, [2, 1, 3]);
            a2 = reshape(a2, rcA(i+2), rx3*rx3);
            a2 = reshape(Ac{i}, rcA(i+1)*rta2, rcA(i+2))*a2;
            a2 = reshape(a2, rcA(i+1), rta2, rx3*rx3);
            a2 = permute(a2, [2,1,3]);
            a2 = reshape(a2, rta2, rcA(i+1)*rx3*rx3);
            ph2 = phfAb{i+1}{L(i+1)+1};
            ph2 = permute(ph2, [1, 3, 2]);
            ph2 = reshape(ph2, n2*n2, rta2);
            a2 = ph2*a2;
            a2 = reshape(a2, n2, n2, rcA(i+1), rx3, rx3);
            a2 = permute(a2, [1, 4, 2, 5, 3]);
            a2 = reshape(a2, n2*rx3, n2*rx3, rcA(i+1));
            
            % All the same for the rhs
            y1 = reshape(lefty1, rx1, rty1, rcy(i+1));
            y1 = permute(y1, [1, 3, 2]);
            y1 = reshape(y1, rx1*rcy(i+1), rty1);
            y1 = y1*(phfyb{i}{L(i)+1}.');
            y1 = reshape(y1, rx1, rcy(i+1), n1);
            y1 = permute(y1, [1, 3, 2]);
            y1 = reshape(y1, rx1*n1, rcy(i+1));
            y2 = phcyr{i+2}.';
            y2 = reshape(yc{i}, rcy(i+1)*rty2, rcy(i+2))*y2;
            y2 = reshape(y2, rcy(i+1), rty2, rx3);
            y2 = permute(y2, [2,1,3]);
            y2 = reshape(y2, rty2, rcy(i+1)*rx3);
            ph2 = phfyb{i+1}{L(i+1)+1};
            y2 = ph2*y2;
            y2 = reshape(y2, n2, rcy(i+1), rx3);
            y2 = permute(y2, [1, 3, 2]);
            y2 = reshape(y2, n2*rx3, rcy(i+1));   
            
            x1 = reshape(xc{i}, rx1*n1, rx2);
            x2 = reshape(xc{i+1}, rx2, n2*rx3).';
            
            if (verb>1)
                fprintf('=core= swp %d, block {%d}, ', swp, i);
            end;
            [u,s,v,r,dx_max,res_max]=local_solve(a1, a2, y1, y2, x1, x2, ...
                rx1, n1, n2, rx3, rcA(i+1), ...
                tol / sqrt(sum(L)) / 2, res_max, dx_max, ...
                local_format, max_full_size, nrestart, gmres_iters, verb);  
            v = v*s;
            
            % kick
            u = reort(u, randn(rx1*n1, kickrank));
            radd = size(u,2)-r;
            v = [v, zeros(n2*rx3, radd)];
            r = size(u,2);
            
            r_max = max(r_max, r);

            % Recompute phis
            % Parts of phis are ready in a1, y1
            a1 = permute(a1, [1, 3, 2]);
            a1 = reshape(a1, rx1*n1*rcA(i+1), rx1*n1);
            curph = a1*u;
            curph = reshape(curph, rx1*n1, rcA(i+1)*r);
            curph = (u')*curph;
            phcAl{i+1} = reshape(curph, r, rcA(i+1), r);
            phcyl{i+1} = (u')*y1;
            
            % Stuff back
            rcx(i+1) = r;
            xc{i} = reshape(u, rx1, n1, r);
            v = reshape(v, n2, rx3, r);
            xc{i+1} = permute(v, [3, 1, 2]);
        end;
    end;
    
    % Exit if we've converged
    if (verb>0)
        fprintf('=rake_solve= swp %d, dx_max: %3.3e, res_max: %3.3e, r_max: %d\n', swp, dx_max, res_max, r_max);
    end;  
    if (res_max<tol)
        break;
    end;
end;

end

function [Phi] = compute_next_Phi(Phi_prev, x, A, y, direction)
% Performs the recurrent Phi (or Psi) matrix computation
% Phi = Phi_prev * (x'Ay).
% If direction is 'lr', computes Psi 
% if direction is 'rl', computes Phi 
% A can be empty, then only x'y is computed.

if (strcmp(direction, 'rl')) 
  % Revert ranks to perform the right-to-left recursion
  x = permute(x, [3, 2, 1]);
  y = permute(y, [3, 2, 1]);
  if (~isempty(A))
    A = permute(A, [4, 2, 3, 1]);
  end
end

rx1 = size(x,1); n = size(x,2); rx2 = size(x,3);
ry1 = size(y,1); ry2 = size(y,3);
if (~isempty(A))
  ra1 = size(A,1); ra2 = size(A,4);
else
  ra1 = 1; ra2 = 1;
end

Phi = reshape(Phi_prev, [rx1*ra1, ry1]);
y = reshape(y, [ry1, n*ry2]);
Phi = Phi*y;	% complexity §\mcommentfont$\mathcal{O}(n  r_x r_A r_y^2)$§
Phi = reshape(Phi, [rx1, ra1, n, ry2]);
Phi = permute(Phi, [2, 3, 1, 4]);
if (~isempty(A))
  Phi = reshape(Phi, [ra1*n, rx1*ry2]);
  A = permute(A, [4, 2, 1, 3]);
  A = reshape(A, [ra2*n, ra1*n]);
  Phi = A*Phi;	% complexity §\mcommentfont$\mathcal{O}(n^2  r_x r_A^2 r_y)$§
  Phi = reshape(Phi, [ra2, n, rx1, ry2]);
end
Phi = permute(Phi, [3, 2, 1, 4]);
Phi = reshape(Phi, [rx1*n, ra2*ry2]);
x = reshape(x, [rx1*n, rx2]);
Phi = (x')*Phi;	% complexity §\mcommentfont$\mathcal{O}(n  r_x^2 r_A r_y)$§
if (~isempty(A))
  Phi = reshape(Phi, [rx2, ra2, ry2]);
end
end

function [A] = core_matrix(core_block, Phi_factor)
% Computes the matrix to use in DMRG by the core, by convolving the core of
% the Tucker block and the last factor Phi matrix over the Tucker rank
% core_block is rc1 x rtuck x rc2, and
% Phi_factor is n1 x rtuck x n2
% A is then rc1 x n1 x n2 x rc2

r1 = size(core_block, 1); rtuck = size(core_block, 2); r2 = size(core_block, 3);
n1 = size(Phi_factor, 1); n2 = size(Phi_factor, 3);

Phi_factor = permute(Phi_factor, [1, 3, 2]);
Phi_factor = reshape(Phi_factor, n1*n2, rtuck); 
A = permute(core_block, [2, 1, 3]);
A = reshape(A, rtuck, r1*r2);
A = Phi_factor*A;
A = reshape(A, n1, n2, r1, r2);
A = permute(A, [3, 1, 2, 4]);

end


function [A] = core_vector(core_block, Phi_factor)
% Computes the vector to use in DMRG by the core, by convolving the core of
% the Tucker block and the last factor Phi matrix over the Tucker rank
% core_block is rc1 x rtuck x rc2, and
% Phi_factor is n1 x rtuck
% A is then rc1 x n1 x rc2

r1 = size(core_block, 1); rtuck = size(core_block, 2); r2 = size(core_block, 3);
n1 = size(Phi_factor, 1);

A = permute(core_block, [2, 1, 3]);
A = reshape(A, rtuck, r1*r2);
A = Phi_factor*A;
A = reshape(A, n1, r1, r2);
A = permute(A, [2, 1, 3]);

end


function [y]=bfun2(B, x, rxm1, m1, m2, rxm3, rxn1, k1, k2, rxn3)
% Computes (B{1} \otimes B{2})x
% B{1} is of sizes rxn1*k1, rxm1*m1, rB
% B{2} is of sizes k2*rxn3, m2*rxm3, rB
rB=size(B{1},3);
x = reshape(x, rxm1*m1, m2*rxm3);
B1 = permute(B{1}, [3 1 2]);
B1 = reshape(B1, rB*rxn1*k1, rxm1*m1);
y = B1*x; % size rB*rxn1*k1,m2*rxm3
y = reshape(y, rB, rxn1*k1, m2*rxm3);
y = permute(y, [3 1 2]);
y = reshape(y, m2*rxm3*rB, rxn1*k1);
B2 = reshape(B{2}, k2*rxn3, m2*rxm3*rB);
y = B2*y; % size k2*rxn3,rxn1*k1
y = reshape(y.', rxn1*k1*k2*rxn3, 1);
end



function [u,s,v,r,dx_max,res_max]=local_solve(a1, a2, y1, y2, x1, x2, ...
    rx1, n1, n2, rx3, ra2, ...
    real_tol, res_max, dx_max, ...
    local_format, max_full_size, nrestart, gmres_iters, verb)

if (strcmp(local_format, 'full'))
    sol_prev = x1*(x2.');
    sol_prev = sol_prev(:);
    rhs = y1*(y2.');
    rhs = rhs(:);
    normy = norm(rhs);
    if (rx1*n1*n2*rx3<max_full_size)
        B = reshape(a1, rx1*n1*rx1*n1, ra2);
        B = B*(reshape(a2, n2*rx3*n2*rx3, ra2).');
        B = reshape(B, rx1*n1,  rx1*n1, n2*rx3, n2*rx3);
        B = permute(B, [1, 3, 2, 4]);
        B = reshape(B, rx1*n1*n2*rx3, rx1*n1*n2*rx3);
        sol = B \ rhs;
        res_true = norm(B*sol-rhs)/normy;
        res_prev = norm(B*sol_prev-rhs)/normy;
    else
        B = cell(2,1);
        B{1} = a1;
        B{2} = a2;
        [sol,flg]=gmres(@(v)bfun2(B, v, rx1, n1, n2, rx3, rx1, n1, n2, rx3), rhs, ...
            nrestart, real_tol, gmres_iters, [], [], sol_prev);
        res_true = norm(bfun2(B, sol, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
        res_prev = norm(bfun2(B, sol_prev, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
    end;
    
    dx = norm(sol-sol_prev)/norm(sol);
    dx_max = max(dx_max, dx);
    res_max = max(res_max, res_prev);
    
    sol = reshape(sol, rx1*n1, n2*rx3);
    [u,s,v]=svd(sol, 'econ');
    s = diag(s);
    r = 1;
    while (r<=numel(s))
        cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
        if (rx1*n1*n2*rx3<max_full_size)
            res = norm(B*cursol(:)-rhs)/normy;
        else
            res = norm(bfun2(B, cursol, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
        end;
        if (res<max(real_tol*2, res_true))
            break;
        end;
        r = r+1;
    end;
    if (verb>1)
        fprintf('dx: %3.3e, res: %3.3e, res_prev: %3.3e, r: %d\n', dx, res, res_prev, r);
    end;
    
    s = diag(s(1:r));
    u = u(:,1:r);
    v = conj(v(:,1:r));
else
    % implement tt-gmres here
end;

end
