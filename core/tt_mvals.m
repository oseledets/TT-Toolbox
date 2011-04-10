function [res] = tt_mvals(a,x,y,nswp)
% [RES] = TT_MVALS(A,X,Y,NSWP)
% This is a function that computes several implicit ALS sweeps for
% a tensor given as a matrix-by-vector product A*X
% Such procedure is supposed to be much cheaper than
% Multiplication and compression afterwards, by a factor of r^2
% A, X are given (in TT format) Y is the initial guess (also in TT format)
% NSWP is the number of ALS iterations (sweeps)
% Output: RES is the improved approximation to A*X
%
% TT Toolbox 1.0, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
res=0;
d=size(x,1);
%prd=tt_mv(a,x); 
%The first step is to compute left-to-right sweep of summation, yielding
%phi-matrices
phi=cell(d,1); %To store phi matrices
nrmf=cell(d,1); %To store their norms
[q,rv]=qr(y{1},0);
y{1}=q;
%A(i1,j1,a) X(j1,a')*Y(i1,a'')
phi{1}=ten_conv(ten_conv(a{1},2,x{1}),1,y{1}); %Is ry x rx x ra  
nrmf{1}=norm(phi{1}(:));
phi{1}=phi{1}./nrmf{1};
%phi{1}=y{1}'*x{1}; %will be ry x rx
for i=2:d-1
    %fprintf('i=%d/%d \n',i,d-1);
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
   %Decipher sizes
   n=size(a{i},1); m=size(a{i},2); ra1=size(a{i},3); ra2=size(a{i},4);
   rx1=size(x{i},2); rx2=size(x{i},3); ry1=size(y{i},2); ry2=size(y{i},3);
   %phi is ry1 x rx1 x ra1, a{i} is n x m x ra1 x ra2 -> n x m x ra2 x rx1
   %keyboard;
   phi{i}=reshape(phi{i-1},[ry1*rx1,ra1])*reshape(permute(a{i},[3,1,2,4]),[ra1,n*m*ra2]);
   %phi{i} is ry1 x rx1 x n x m x ra2, over m,rx1
   phi{i}=reshape(x{i},[m*rx1,rx2])'*reshape(permute(reshape(phi{i},[ry1,rx1,n,m,ra2]),[4,2,3,1,5]),[m*rx1,n*ry1*ra2]); %m x rx1 x n x ry1 x ra2
   %phi{i} is rx2 x n x ry1 x ra2
%   if ( i == 190 )
%       keyboard
%   end
   phi{i}=reshape(y{i},[n*ry1,ry2])'*reshape(permute(reshape(phi{i},[rx2,n*ry1,ra2]),[2,1,3]),[n*ry1,rx2*ra2]);
   phi{i}=reshape(phi{i},[ry2,rx2,ra2]);
   nrmf{i}=norm(phi{i}(:));
   phi{i}=phi{i}./nrmf{i};
end
%Now we are able to recompute Y_d
%fprintf('Last phi-norm %e \n',norm(phi{d-1}(:)));
for s=1:nswp
%Y(id,a'')=\sum_{jd,a,a'} A(id,jd,a)*X(jd,a')*phi(a,a',a'') 
n=size(a{d},1); m=size(a{d},2); ra=size(a{d},3);
rx=size(x{d},2); ry=size(phi{d-1},1);
y{d}=reshape(reshape(permute(a{d},[1,3,2]),[n*ra,m])*x{d},[n,ra*rx])*reshape(permute(phi{d-1},[1,3,2]),[ry,ra*rx])';
%fprintf('y-mv \n');
%y{d}

%res=y;
%return;
%Now we can start right-to-left sweep that combines computation of 
%phi matrices with Y right-to-left orthogonalization
[q,rv]=qr(y{d},0);
y{d}=q;
ry=size(y{d},2);
%a{d} is n x m x ra, x{d} is m x rx, y{d} is m x ry, phi{d} will be ry x rx x ra
%phi(a,a',a'') = A(i,j,a'')*X(j,a')*Y(i,a)
phi{d}=y{d}'*reshape(reshape(permute(a{d},[1,3,2]),[n*ra,m])*x{d},[n,ra*rx]); %is n x ra x rx
phi{d}=permute(reshape(phi{d},[ry,ra,rx]),[1,3,2]);
nrmf{d}=norm(phi{d}(:));
phi{d}=phi{d}./nrmf{d};
%fprintf('mv: \n');
%phi{d}
%res=y;
%return
%fprintf('mv phi{%d} norm: %e \n',d,norm(phi{d}(:)));
for i=d-1:-1:2
    %fprintf('i=%d/%d \n',i,d-1);
   %First solve for y{i} then orthogonalize it
   n=size(a{i},1); m=size(a{i},2); ra1=size(a{i},3); ra2=size(a{i},4);
   rx1=size(x{i},2); rx2=size(x{i},3); ry1=size(y{i},2); ry2=size(phi{i+1},1);
   y{i}=reshape(a{i},[n*m*ra1,ra2])*reshape(phi{i+1},[ry2*rx2,ra2])';
   y{i}=reshape(y{i},[n,m,ra1,ry2,rx2]);
   %is n x m x ra1 x ry2 x rx2
   phi0=reshape(permute(x{i},[1,3,2]),[m*rx2,rx1])'*reshape(permute(y{i},[2,5,3,1,4]),[m*rx2,ra1*n*ry2]);
   %is rx1 x ra1 x n x ry2 phi0 is a convolution of (A*X) with phi{i+1}
   %over (ra2,rx2), third index, in "old" life it was n x R x ry2, R = ra2
   %x rx2
   phi0=reshape(phi0,[rx1*ra1,n*ry2]);
   %fprintf('mv \n'); permute(reshape(phi0',[n,ry2,rx1*ra1]),[1,3,2])
   %   return
   
   y{i}=reshape(phi{i-1},[ry1,rx1*ra1])*phi0;
   
   %is ry1 x n x ry2
   y{i}=permute(reshape(y{i},[ry1,n,ry2]),[2,1,3]);
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
    phi{i}=permute(reshape(phi0*reshape(permute(y{i},[1,3,2]),[n*ry2,ry1]),[rx1,ra1,ry1]),[3,1,2]); % rx1 x ra1 x ry1 
   nrmf{i}=norm(phi{i}(:));
   phi{i}=phi{i}./nrmf{i};
    %fprintf('mv phi{%i} norm: %e \n',i,norm(phi{i}(:)));
end

%fprintf('mv: \n');
%for i=2:d
%  fprintf('%e \n',norm(phi{i}(:)));
%end
%y{1} = x{1}*phi{2}';
%This one was utterly forgotten
%A(i,j,ra)*X(j,rx)*phi(ry,rx,ra)
n=size(a{1},1); m=size(a{1},2); ra=size(a{1},3);
ry=size(phi{2},1);rx=size(phi{2},2);
%keyboard;
y{1}=reshape(reshape(permute(a{1},[1,3,2]),[n*ra,m])*x{1},[n,rx*ra])*...
reshape(permute(phi{2},[1,3,2]),[ry,rx*ra])'; %is n x ra x rx
%then n x ry
%res=y;
%phi{2}
%return;
%Now (finally!) we have this left-to-right sweep
[q,rv]=qr(y{1},0);
y{1}=q;
%X(i1,a)*Y(i1,b)  
%A(n,m,ra)*X(m,rx)*Y(n,ry)
phi{1}=ten_conv(ten_conv(a{1},1,y{1}),2,x{1});
nrmf{1}=norm(phi{1}(:));
phi{1}=phi{1}./nrmf{1};
%phi{1}=permute(phi{1},[1,3,2]);
%phi{1}
for i=2:d-1
   %fprintf('i=%d/%d \n',i,d-1);
   %A(n,m,ra1,ra2)*X(m,rx1,rx2)*phi(ry1,rx1,ra1) -> 
   n=size(a{i},1); m=size(a{i},2); ra1=size(a{i},3); ra2=size(a{i},4);
   rx1=size(x{i},2); rx2=size(x{i},3); ry1=size(phi{i-1},1);
   phi0=reshape(permute(x{i},[2,1,3]),[rx1,m*rx2])'*reshape(permute(phi{i-1},[2,1,3]),[rx1,ry1*ra1]);
   %phi0 is m x rx2 x ry1 x ra1, convolve with A(n,m,ra1,ra2)
   phi0=reshape(permute(reshape(phi0,[m,rx2,ry1,ra1]),[1,4,2,3]),[m*ra1,rx2*ry1]);
   phi0=reshape(permute(a{i},[2,3,1,4]),[m*ra1,n*ra2])'*phi0;
   %phi0 is n x ra2 x rx2 x ry1
   %phi{i+1}(ry2,rx2,ra2)
   %r^2 x mr x mr x nr -> nmr^4
   ry2=size(phi{i+1},1);
   y{i}=reshape(permute(reshape(phi0,[n,ra2,rx2,ry1]),[1,4,3,2]),[n*ry1,rx2*ra2])*...
       reshape(phi{i+1},[ry2,rx2*ra2])';
   y{i}=reshape(y{i},[n,ry1,ry2]);
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
   phi{i}=reshape(y{i},[n*ry1,ry2])'* ...
       reshape(permute(reshape(phi0,[n,ra2,rx2,ry1]),[1,4,3,2]),[n*ry1,rx2*ra2]);
   phi{i}=reshape(phi{i},[ry2,rx2,ra2]);
   nrmf{i}=norm(phi{i}(:));
   phi{i}=phi{i}./nrmf{i};
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
rx=size(x{d},2); ry=size(phi{d-1},1);
y{d}=reshape(reshape(permute(a{d},[1,3,2]),[n*ra,m])*x{d},[n,ra*rx])*reshape(permute(phi{d-1},[1,3,2]),[ry,ra*rx])';
%wi=reshape(permute(phi{d-1},[1,3,2]),[ry,ra*rx]);
%y{d}=prd{d}*wi';
%y{d}= x{d}*phi{d-1}';
for i=2:d
   y{i}=y{i}.*nrmf{i-1};
end
res=y;
return;
end
