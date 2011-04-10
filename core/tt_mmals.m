function [res] = tt_mmals(a,x,y,nswp)
% [RES] = TT_MMALS(A,X,Y,NSWP)
% This is a function that computes several implicit ALS sweeps for
% a tensor given as a matrix-by-matrix product A*X
% Such procedure is supposed to be much cheaper than
% Multiplication and compression afterwards, by a factor of r^2
% A, X are given (in TT format) Y is the initial guess (also in TT format)
% NSWP is the number of ALS iterations (sweeps)
% Output: RES is the improved approximation to A*X
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
res=0;
d=size(x,1);
%The first step is to compute left-to-right sweep of summation, yielding
%phi-matrices
phi=cell(d,1); %To store phi matrices
nrmf=cell(d,1); %To store their norms
n=size(y{1},1); s=size(y{1},2); ry=size(y{1},3);
[q,rv]=qr(reshape(y{1},[n*s,ry]),0);
y{1}=q; ry=size(q,2); y{1}=reshape(y{1},[n,s,ry]);
%A(i1,j1,a) X(j1,k1,a')*Y(i1,k1,a'')
n=size(a{1},1); m=size(a{1},2); ra=size(a{1},3); s=size(x{1},2);
rx=size(x{1},3); ry=size(y{1},3);
phi{1}=ten_conv(a{1},2,reshape(x{1},[m,s*rx])); %is n x (s) x rx x ra
phi{1}=reshape(y{1},[n*s,ry])'*reshape(phi{1},[n*s,rx*ra]);
%phi{1} is ry x rx x ra
phi{1}=reshape(phi{1},[ry,rx,ra]); %Is ry x rx x ra  
nrmf{1}=norm(phi{1}(:));
%phi{1}=phi{1}./nrmf{1}; %Disabled for testing purposes
for i=2:d-1
    %fprintf('i=%d/%d \n',i,d-1);
    n=size(y{i},1); s=size(y{i},2); ry1=size(y{i},3); ry2=size(y{i},4);
    y{i}=reshape(y{i},[n*s,ry1,ry2]);
    y{i} = ten_conv(y{i},2,rv');
    ncur=size(y{i},1);
    r2=size(y{i},2);
    r3=size(y{i},3);
    core=permute(y{i},[1,2,3]);
    core=reshape(core,[ncur*r2,r3]);
    [y{i},rv]=qr(core,0);
    rv=rv./norm(rv(:));
    rnew=min(r3,ncur*r2);
    y{i}=reshape(y{i},[ncur,r2,rnew]);
    y{i}=permute(y{i},[1,2,3]);
    y{i}=reshape(y{i},[n,s,r2,rnew]);
    %Decipher sizes
   n=size(a{i},1); m=size(a{i},2); s=size(x{i},2); 
   ra1=size(a{i},3); ra2=size(a{i},4);
   rx1=size(x{i},3); rx2=size(x{i},4); ry1=size(y{i},3); ry2=size(y{i},4);
   %phi is ry1 x rx1 x ra1, a{i} is n x m x ra1 x ra2 -> n x m x ra2 x rx1
   %keyboard;
   phi{i}=reshape(phi{i-1},[ry1*rx1,ra1])*reshape(permute(a{i},[3,1,2,4]),[ra1,n*m*ra2]);
   %phi{i} is ry1 x rx1 x n x m x ra2, over m,rx1
   %x{i} is m x s x rx1 x rx2
   %After reshape phi{i} is m x rx1 x n*ry1*ra2
   
   phi{i}=reshape(permute(x{i},[1,3,2,4]),[m*rx1,s*rx2])'*reshape(permute(reshape(phi{i},[ry1,rx1,n,m,ra2]),[4,2,3,1,5]),[m*rx1,n*ry1*ra2]); %m x rx1 x n x ry1 x ra2
   %phi{i} is s x rx2 x n x ry1 x ra2,
   %convolve over s,n,ry1
   phi{i}=reshape(phi{i},[s,rx2,n,ry1,ra2]);
   phi{i}=permute(phi{i},[3,1,4,2,5]); %is n x s x ry1 x rx2 x ra2
   phi{i}=reshape(y{i},[n*s*ry1,ry2])'*reshape(phi{i},[n*s*ry1,rx2*ra2]);
   %phi{i}=reshape(phi{i},
%   if ( i == 190 )
%       keyboard
%   end
   %phi{i}=reshape(y{i},[n*ry1,ry2])'*reshape(permute(reshape(phi{i},[rx2,n*ry1,ra2]),[2,1,3]),[n*ry1,rx2*ra2]);
   phi{i}=reshape(phi{i},[ry2,rx2,ra2]);
   nrmf{i}=norm(phi{i}(:));
   %phi{i}=phi{i}./nrmf{i}; Disabled for testing purposes
end
%Now we are able to recompute Y_d
%fprintf('Last phi-norm %e \n',norm(phi{d-1}(:)));
for s0=1:nswp
%Y(id,kd,a'')=\sum_{jd,a,a'} A(id,jd,a)*X(jd,kd,a')*phi(a,a',a'') 
n=size(a{d},1); m=size(a{d},2); ra=size(a{d},3);
rx=size(x{d},3); s=size(x{d},2); ry=size(phi{d-1},1);
%n x m x ra m x s x rx
y{d}=reshape(permute(a{d},[1,3,2]),[n*ra,m])*reshape(x{d},[m,s*rx]);
%y{d} is n x ra x s x rx
y{d}=reshape(permute(reshape(y{d},[n,ra,s,rx]),[1,3,4,2]),[n*s,rx*ra])*...
 reshape(phi{d-1},[ry,rx*ra])';
y{d}=reshape(y{d},[n,s,ry]);
%keyboard;
%y{d}=reshape(reshape(permute(a{d},[1,3,2]),[n*ra,m])*x{d},[n,ra*rx])*...
%reshape(permute(phi{d-1},[1,3,2]),[ry,ra*rx])';
%fprintf('y-mv \n');
%y{d}

%res=y;
%return;
%Now we can start right-to-left sweep that combines computation of 
%phi matrices with Y right-to-left orthogonalization
n=size(y{d},1); s=size(y{d},2); ry=size(y{d},3);
y{d}=reshape(y{d},[n*s,ry]);
[q,rv]=qr(y{d},0);
y{d}=q;
ry=size(y{d},2);
y{d}=reshape(y{d},[n,s,ry]);
%a{d} is n x m x ra, x{d} is m x rx, y{d} is m x ry, phi{d} will be ry x rx x ra
%phi(a,a',a'') = A(i,j,a'')*X(j,k,a')*Y(i,k,a)

rx=size(x{d},3); m=size(x{d},1); s=size(y{d},2);
phi{d}=reshape(permute(a{d},[1,3,2]),[n*ra,m])*reshape(x{d},[m,s*rx]);
%phi{d} is n x ra x s x rx, (n,s)
phi{d}=reshape(y{d},[n*s,ry])'*reshape(permute(reshape(phi{d},...
    [n,ra,s,rx]),[1,3,4,2]),[n*s,rx*ra]);
    
phi{d}=reshape(phi{d},[ry,rx,ra]);
%phi{d}=y{d}'*reshape(reshape(permute(a{d},[1,3,2]),[n*ra,m])*x{d},[n,ra*rx]); %is n x ra x rx
%phi{d}=permute(reshape(phi{d},[ry,ra,rx]),[1,3,2]);


%nrmf{d}=norm(phi{d}(:));
%phi{d}=phi{d}./nrmf{d};


%fprintf('mv: \n');
%phi{d}
%res=y;
%return
%fprintf('mv phi{%d} norm: %e \n',d,norm(phi{d}(:)));
for i=d-1:-1:2
    %fprintf('i=%d/%d \n',i,d-1);
   %First solve for y{i} then orthogonalize it
   n=size(a{i},1); m=size(a{i},2); ra1=size(a{i},3); ra2=size(a{i},4);
   s=size(x{i},2); 
   rx1=size(x{i},3); rx2=size(x{i},4); ry1=size(y{i},3); ry2=size(phi{i+1},1);
   y{i}=reshape(a{i},[n*m*ra1,ra2])*reshape(phi{i+1},[ry2*rx2,ra2])';
   y{i}=reshape(y{i},[n,m,ra1,ry2,rx2]);
   %is n x m x ra1 x ry2 x rx2
   %Convolve over m,rx2,
   %phi0=reshape(permute(x{i},[1,3,2]),[m*rx2,rx1])'*reshape(permute(y{i},[2,5,3,1,4]),[m*rx2,ra1*n*ry2]);
   phi0=reshape(permute(y{i},[1,4,3,5,2]),[n*ry2*ra1,rx2*m])*reshape(permute(x{i},[4,1,2,3]),[rx2*m,s*rx1]);
   %phi0 is n x ry2 x ra1 x s x rx1 
   phi0=reshape(phi0,[n,ry2,ra1,s,rx1]);
   phi0=permute(phi0,[5,3,1,4,2]);
   %phi0=reshape(phi0,[n*s*ry2,rx1*ra1]);
   %y{i}
   %is rx1 x ra1 x n x ry2 phi0 is a convolution of (A*X) with phi{i+1}
   %over (ra2,rx2), third index, in "old" life it was n x R x ry2, R = ra2
   %x rx2
   phi0=reshape(phi0,[rx1*ra1,n*s*ry2]);
   %fprintf('mv \n'); permute(reshape(phi0',[n,ry2,rx1*ra1]),[1,3,2])
   %   return
   
   y{i}=reshape(phi{i-1},[ry1,rx1*ra1])*phi0;
   
   %is ry1 x n x ry2
   ncur=n*s;
   y{i}=permute(reshape(y{i},[ry1,ncur,ry2]),[2,1,3]);
   %wi=reshape(permute(phi{i+1},[1,3,2]),[ry2,ra2*rx2]);
   %qi=reshape(permute(phi{i-1},[1,3,2]),[ry1,ra1*rx1]);
   %y{i}=ten_conv(prd{i},3,wi');
   %y{i}=ten_conv(y{i},2,qi');
   %Orthogonalize y
    ncur=size(y{i},1);
    r2=size(y{i},2);
    r3=size(y{i},3);
    core=permute(y{i},[1,3,2]);
    core=reshape(core,[ncur*r3,r2]);
    [y{i},rv]=qr(core,0);
    rnew=min(r2,ncur*r3);
    y{i}=reshape(y{i},[ncur,r3,rnew]);
    y{i}=permute(y{i},[1,3,2]);
  % fprintf('mv: \n'); y{i}
  % return;
  
    %Compute new phi{i}, which is what?
   %it was convolution of old X with phi{i+1} from the right (we did it before!!)
   %And convolution of it with new y{i}
  
   % keyboard;
    %phi0=reshape(x{i},[ncur*rx1,rx2])*phi{i+1}';
    %phi0 is rx1 x ra1 x n x ry2, y{i} is n x ry1 x ry2
    %phi{i} is ry1 x rx1 x ra1, convolve over n & ry2
    ry1=size(y{i},2); ry2=size(y{i},3);    
    phi{i}=permute(reshape(phi0*reshape(permute(y{i},[1,3,2]),[ncur*ry2,ry1]),...
    [rx1,ra1,ry1]),[3,1,2]); % rx1 x ra1 x ry1 
   nrmf{i}=norm(phi{i}(:));
   %phi{i}=phi{i}./nrmf{i};
    %fprintf('mv phi{%i} norm: %e \n',i,norm(phi{i}(:)));
   y{i}=reshape(y{i},[n,s,ry1,ry2]);
end
%keyboard;
%fprintf('mv: \n');
%for i=2:d
%  fprintf('%e \n',norm(phi{i}(:)));
%end
%y{1} = x{1}*phi{2}';
%This one was utterly forgotten
%A(n,m,ra)*X(m,s,rx)*phi(ry,rx,ra)
n=size(a{1},1); m=size(a{1},2); ra=size(a{1},3);
s=size(x{1},2);
ry=size(phi{2},1);rx=size(phi{2},2);
%keyboard;
y{1}=reshape(permute(a{1},[1,3,2]),[n*ra,m])*reshape(x{1},[m,s*rx]);%,[n*s,rx*ra])*...
%y{1} is n x ra x s x rx
y{1}=reshape(permute(reshape(y{1},[n,ra,s,rx]),[1,3,4,2]),[n*s,rx*ra])*...
reshape(phi{2},[ry,rx*ra])'; %is n x ra x rx
y{1}=reshape(y{1},[n,s,ry]);
%then n x ry
%res=y;
%phi{2}
%return;
%Now (finally!) we have this left-to-right sweep
y{1}=reshape(y{1},[n*s,ry]);
[q,rv]=qr(y{1},0);
y{1}=q;
y{1}=reshape(y{1},[n,s,size(y{1},2)]);
%%X(i1,a)*Y(i1,b)  
%A(n,m,ra)*X(m,s,rx)*Y(n,s,ry)
n=size(a{1},1); m=size(a{1},2); ra=size(a{1},3);
s=size(x{1},2); rx=size(x{1},3); ry=size(y{1},3);
phi{1}=reshape(permute(a{1},[1,3,2]),[n*ra,m])*reshape(x{1},[m,s*rx]);
%phi{1} is n x ra x s x rx
phi{1}=reshape(permute(reshape(phi{1},[n,ra,s,rx]),[1,3,4,2]),[n*s,rx*ra])'*...
    reshape(y{1},[n*s,ry]);
phi{1}=permute(reshape(phi{1},[rx,ra,ry]),[3,1,2]);
%phi{1}=ten_conv(ten_conv(a{1},1,y{1}),2,x{1});
%nrmf{1}=norm(phi{1}(:)); Fix it later
%phi{1}=phi{1}./nrmf{1};
%phi{1}=permute(phi{1},[1,3,2]);
%phi{1}
for i=2:d-1
   %fprintf('i=%d/%d \n',i,d-1);
   %A(n,m,ra1,ra2)*X(m,s,rx1,rx2)*phi(ry1,rx1,ra1) -> 
   n=size(a{i},1); m=size(a{i},2); ra1=size(a{i},3); ra2=size(a{i},4);
   s=size(x{i},2);
   rx1=size(x{i},3); rx2=size(x{i},4); ry1=size(phi{i-1},1);
   phi0=reshape(permute(x{i},[3,1,2,4]),[rx1,m*s*rx2])'*reshape(permute(phi{i-1},[2,1,3]),[rx1,ry1*ra1]);
   %phi0 is m x s x rx2 x ry1 x ra1, convolve with A(n,m,ra1,ra2)
   phi0=reshape(permute(reshape(phi0,[m,s,rx2,ry1,ra1]),[1,5,2,3,4]),[m*ra1,s*rx2*ry1]);
   phi0=reshape(permute(a{i},[2,3,1,4]),[m*ra1,n*ra2])'*phi0;
   phi0=permute(reshape(phi0,[n,ra2,s,rx2,ry1]),[1,3,2,4,5]);
   %phi0 is n x ra2 x s x rx2 x ry1
   %phi{i+1}(ry2,rx2,ra2)
   %r^2 x mr x mr x nr -> nmr^4
   ncur=n*s;
   ry2=size(phi{i+1},1);
   y{i}=reshape(permute(reshape(phi0,[ncur,ra2,rx2,ry1]),[1,4,3,2]),[ncur*ry1,rx2*ra2])*...
       reshape(phi{i+1},[ry2,rx2*ra2])';
   y{i}=reshape(y{i},[ncur,ry1,ry2]);
   %norm(y{i}(:))
   %return
   %y{i}=ten_conv(x{i},3,phi{i+1}');
   %y{i}=ten_conv(y{i},2,phi{i-1}');
   %Orthogonalize new core
   ncur=size(y{i},1);
   r2=size(y{i},2);
   r3=size(y{i},3);
   core=permute(y{i},[1,2,3]);
   core=reshape(core,[ncur*r2,r3]);
   [y{i},rv]=qr(core,0);
   rnew=min(r3,ncur*r2);
   y{i}=reshape(y{i},[ncur,r2,rnew]);
   y{i}=permute(y{i},[1,2,3]);
   %Compute new phi
   %n=size(x{i},1); rx1=size(x{i},2);  rx2=size(x{i},3);
   ry1=size(y{i},2); ry2=size(y{i},3);
   %phi{i}=phi{i-1}*reshape(permute(x{i}, [2,1,3]),[rx1,n*rx2]);
   %How to: convolve x over phi{i-1}, then convolve the result with new y.
   %phi0 is n x ra2 x rx2 x ry1, y is n x ry1 x ry2
   %phi{i} is ry2 x rx2 x ra2
   phi{i}=reshape(y{i},[ncur*ry1,ry2])'* ...
       reshape(permute(reshape(phi0,[ncur,ra2,rx2,ry1]),[1,4,3,2]),[ncur*ry1,rx2*ra2]);
   phi{i}=reshape(phi{i},[ry2,rx2,ra2]);
   %nrmf{i}=norm(phi{i}(:));
   %phi{i}=phi{i}./nrmf{i};
   y{i}=reshape(y{i},[n,s,ry1,ry2]);
   
   %pri=reshape(permute(phi{i-1},[1,3,2]),[ry1,ra1*rx1]);
   %n=size(prd{i},1); rx1=size(prd{i},2); rx2=size(prd{i},3);
   %phi{i}=pri*reshape(permute(prd{i}, [2,1,3]),[rx1,n*rx2]);
   %phi is now ry1 x n x rx2, core_y is n x ry1 x ry2  
   %phi{i}=permute(reshape(phi{i},[ry1,n,rx2]),[2,1,3]); % n x ry1 x rx2
   %phi{i}=reshape(y{i},[n*ry1,ry2])'*reshape(phi{i},[n*ry1,rx2]);
   %rx1=size(x{i},2); rx2=size(x{i},3);
   %phi{i}=permute(reshape(phi{i},[ry2,ra2,rx2]),[1,3,2]);
end
%Now we are able to recompute Y_d
%mat_x{d} is X(i_d,rx) and phi is  ry x rx
end
n=size(a{d},1); m=size(a{d},2); ra=size(a{d},3);
rx=size(x{d},3); s=size(x{d},2); ry=size(phi{d-1},1);
%n x m x ra m x s x rx
y{d}=reshape(permute(a{d},[1,3,2]),[n*ra,m])*reshape(x{d},[m,s*rx]);
%y{d} is n x ra x s x rx
y{d}=reshape(permute(reshape(y{d},[n,ra,s,rx]),[1,3,4,2]),[n*s,rx*ra])*...
 reshape(phi{d-1},[ry,rx*ra])';
y{d}=reshape(y{d},[n,s,ry]);%wi=reshape(permute(phi{d-1},[1,3,2]),[ry,ra*rx]);
%y{d}=prd{d}*wi';
%y{d}= x{d}*phi{d-1}';
%for i=2:d
%   y{i}=y{i}.*nrmf{i-1};
%end
res=y;
return;
end
