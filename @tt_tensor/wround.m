function [tt] = wround(tt, w, eps, x0, rmax)
%[TT]=WROUND(TT,W,[EPS],[X0],[RMAX])
%Approximate TT-tensor with relative accuracy EPS in the norm |||x||| = ||Wx||
%With initial approximation x0 and maximal rank RMAX
%W can be either a symmetric tt_matrix or a tt_tensor. In the latter case it
%is converted to a diagonal matrix (to imply Hadamard product with W).
%eps = 1e-14;
if ( nargin == 2 || isempty(eps) )
   eps=1e-14;
end
if ( nargin == 3 || isempty(x0) )
  x0=tt;
end
if ( nargin == 4 || isempty(rmax) )
  rmax=prod(size(tt));
end
tt=x0;

%Parameters section
nswp = 1000;
if (~isa(w, 'tt_matrix'))
    w = diag(w);
end;

% MATLAB's object-oriented system sucks!!
% We have to do the following crap:
w_vec = w.tt;
rw = w_vec.r;
w_pos = w_vec.ps;
w_cr = w_vec.core;


swp = 1;
while ( swp <= nswp ) 
d=tt.d;
n=tt.n;
r=tt.r;
pos=tt.ps;
cr=tt.core;

pos1=1;
nrm=zeros(d,1);
core0=cr(1:r(1)*n(1)*r(2));
phi = cell(d,1);
phi{1}=1; % the first left phi
%Orthogonalization from left-to-right and compute phi=X^T diag(W) X
for i=1:d-1
   core0=reshape(core0,[r(i)*n(i),r(i+1)]);
   [core0,ru]=qr(core0,0); nrm(i+1)=norm(ru,'fro');
   ru=ru./nrm(i+1);
   core1=cr(pos(i+1):pos(i+2)-1);
   core1=reshape(core1,[r(i+1),n(i+1)*r(i+2)]);
   core1=ru*core1;
   r(i+1)=size(core0,2);
   x1 = reshape(core0, r(i), n(i)*r(i+1)); % we'll need it for phi
   cr(pos1:pos1-1+r(i)*n(i)*r(i+1))=core0(:);
   cr(pos1+r(i)*n(i)*r(i+1):pos1+r(i)*n(i)*r(i+1)+r(i+1)*n(i+1)*r(i+2)-1)=core1(:);
   core0=core1;
   pos1=pos1+r(i)*n(i)*r(i+1);
   
   % now, update phi
   phi{i+1} = reshape(phi{i}, r(i)*rw(i), r(i));
   phi{i+1}=phi{i+1}*x1; % size rx1*rw1, n1*rx2
   phi{i+1}=reshape(phi{i+1}, r(i), rw(i), n(i), r(i+1));
   phi{i+1}=permute(phi{i+1}, [1 4 2 3]);
   phi{i+1}=reshape(phi{i+1}, r(i)*r(i+1), rw(i)*n(i));
   wm1 = w_cr(w_pos(i):w_pos(i+1)-1);
   wm1 = reshape(wm1, rw(i), n(i), n(i), rw(i+1));
   wm1 = reshape(permute(wm1, [1 2 4 3]), rw(i)*n(i), rw(i+1)*n(i));
   phi{i+1}=phi{i+1}*wm1; % size rx1*rx2, rw2*n
   phi{i+1}=reshape(phi{i+1}, r(i), r(i+1), rw(i+1), n(i));
   phi{i+1}=reshape(permute(phi{i+1}, [2 3 1 4]), r(i+1)*rw(i+1), r(i)*n(i));
   phi{i+1}=phi{i+1}*reshape(x1, r(i)*n(i), r(i+1)); % size r(i+1)*rw(i+1),r(i+1)
   phi{i+1}=reshape(phi{i+1}, r(i+1), rw(i+1), r(i+1));
end
pos1=pos1+r(d)*n(d)*r(d+1)-1;
cr=cr(1:pos1); %Truncate storage if required
ep=eps/sqrt(d-1);
pos=cumsum([1;n.*r(1:d).*r(2:d+1)]); 

phi{d}=1; % the last right phi
for i=d:-1:2        
    core2 = cr(pos1-r(i)*n(i)*r(i+1)+1:pos1);
    core1 = cr(pos(i-1):pos(i)-1);
    
    % supercore
    score = reshape(core1, r(i-1)*n(i-1), r(i))*reshape(core2, r(i), n(i)*r(i+1));
    % matrix
    wm1 = w_cr(w_pos(i-1):w_pos(i)-1);
    wm1 = reshape(wm1, rw(i-1), n(i-1), n(i-1), rw(i));
    wm1 = permute(wm1, [1,4,2,3]);
    wm2 = w_cr(w_pos(i):w_pos(i+1)-1);
    wm2 = reshape(wm2, rw(i), n(i), n(i), rw(i+1));
    wm2 = permute(wm2, [4 1 2 3]);
    mcore = cell(2,1);
    % first mblock
    mcore{1} = reshape(permute(phi{i-1}, [1 3 2]), r(i-1)*r(i-1), rw(i-1));
    mcore{1} = mcore{1}*reshape(wm1, rw(i-1), rw(i)*n(i-1)*n(i-1)); % size rx1*rx1, rw2*n1*m1
    mcore{1} = reshape(mcore{1}, r(i-1), r(i-1), rw(i), n(i-1),n(i-1));
    mcore{1} = permute(mcore{1}, [1 4 2 5 3]);
    mcore{1} = reshape(mcore{1}, r(i-1)*n(i-1), r(i-1)*n(i-1), rw(i));
    % second mblock
    mcore{2} = reshape(permute(phi{i}, [1 3 2]), r(i+1)*r(i+1), rw(i+1));
    mcore{2} = mcore{2}*reshape(wm2, rw(i+1), rw(i)*n(i)*n(i)); % size rx3*rx3, rw2*n2*m2
    mcore{2} = reshape(mcore{2}, r(i+1), r(i+1), rw(i), n(i),n(i));
    mcore{2} = permute(mcore{2}, [4 1 5 2 3]);
    mcore{2} = reshape(mcore{2}, n(i)*r(i+1), n(i)*r(i+1), rw(i));
    
    rhs = tt_mat_full_vec(mcore, score(:));
    [u,s,v]=svd(score, 'econ');
    for rnew=1:min(size(u,2), rmax)
        curx = u(:,1:rnew)*s(1:rnew,1:rnew)*(v(:,1:rnew)');
        resid = norm(tt_mat_full_vec(mcore, curx(:))-rhs)/norm(rhs);
        err = norm(curx(:)-score(:))/norm(score(:));
        fprintf('block %d-%d, rank=%d, resid=%3.3e, err=%3.3e\n', i, i-1, rnew, resid, err);
        if (resid<ep)        
            break;
        end;
    end;
%     % bin search for the rank - it sucks!
%     r0 = 1; r1 = min(size(u,2), rmax); rnew = 2;
%     while (rnew~=r0)&&(rnew~=r1)
%         rnew = floor((r1+r0)*0.5);
%         curx = u(:,1:rnew)*s(1:rnew,1:rnew)*(v(:,1:rnew)');
%         resid = norm(tt_mat_full_vec(mcore, curx(:))-rhs)/norm(rhs);
%         fprintf('block %d-%d, rank=%d, resid=%3.3e\n', i, i-1, rnew, resid);
%         if (resid>ep)
%             r0 = rnew+1;
%         else
%             r1 = rnew;
%         end;
%     end;
    r(i)=rnew;
    core1 = u(:,1:rnew)*s(1:rnew,1:rnew);
    core2 = v(:,1:rnew)';
    x2 = core2; % for phi
    
    cr(pos1-r(i)*n(i)*r(i+1)+1:pos1)=core2(:);
    cr(pos1-r(i)*n(i)*r(i+1)-r(i-1)*n(i-1)*r(i)+1:pos1-r(i)*n(i)*r(i+1))=core1(:);
    pos1=pos1-r(i)*n(i)*r(i+1);
    
    % update phi
    phi{i-1} = reshape(phi{i}, r(i+1), rw(i+1)*r(i+1));
    phi{i-1}=reshape(x2, r(i)*n(i), r(i+1))*phi{i-1}; % size r2*n2, rw3*r3
    phi{i-1}=reshape(phi{i-1}, r(i), n(i), rw(i+1), r(i+1));
    phi{i-1}=permute(phi{i-1}, [3 2 1 4]);
    phi{i-1}=reshape(phi{i-1}, rw(i+1)*n(i), r(i)*r(i+1));    
    wm2 = permute(wm2, [1 3 2 4]);
    wm2 = reshape(wm2, rw(i+1)*n(i), rw(i)*n(i));
    phi{i-1}=(wm2.')*phi{i-1}; % size rw2*n2, r2*r3
    phi{i-1}=reshape(phi{i-1}, rw(i), n(i), r(i), r(i+1));
    phi{i-1}=permute(phi{i-1}, [2 4 1 3]);
    phi{i-1}=reshape(phi{i-1}, n(i)*r(i+1), rw(i)*r(i));
    phi{i-1}=reshape(x2, r(i), n(i)*r(i+1))*phi{i-1};
    phi{i-1}=reshape(phi{i-1}, r(i), rw(i), r(i));
end;

 pos1=pos1-r(1)*n(1)*r(2);
 cr=cr(pos1+1:numel(cr)); %Truncate unwanted elements;
 tt.r=r;
 tt.ps=cumsum([1;tt.n.*tt.r(1:d).*tt.r(2:d+1)]);
 pp=cr(1:r(1)*n(1)*r(2));
 nrm(1)=norm(pp,'fro');
 cr(1:r(1)*n(1)*r(2))=pp./nrm(1);
 %Now a simple trick: balance the product of numbers;
 %All cores are orthogonal except the first one. Thus, we know the norm
 nrm0=sum(log(abs(nrm))); 
 nrm0=nrm0/d; nrm0=exp(nrm0);
 %Construct normalization of norm
 for i=1:d-1
   nrm(i+1)=nrm(i+1)*nrm(i)/nrm0;
   nrm(i)=nrm0;
 end
 %Finally redistribute the norm
 ps=tt.ps;
 for i=1:d
    core1=cr(ps(i):ps(i+1)-1);
    core1=core1*nrm(i);
    cr(ps(i):ps(i+1)-1)=core1;
 end
 tt.core=cr;
 swp=swp+1;
end
end