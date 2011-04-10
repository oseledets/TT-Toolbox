function [res] = tt_als(x,y,nswp)
% [RES]=TT_ALS(X,Y,NSWP)
% This is a special function that computes several ALS sweeps for
% TT format. 
% X is given (in TT format), Y is the initial guess (also in TT format)
% Output: RES --- Improved TT approximation to X
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

[q,rv]=qr(y{1},0);
y{1}=q;
%X(i1,a)*Y(i1,b)  
phi{1}=y{1}'*x{1}; %will be ry x rx
for i=2:d-1
    %fprintf('i=%d/%d \n',i,d-1);
    y{i} = ten_conv(y{i},2,rv');
    ncur=size(y{i},1);
    r2=size(y{i},2);
    r3=size(y{i},3);
    core=permute(y{i},[1,2,3]);
    core=reshape(core,[ncur*r2,r3]);
    [y{i},rv]=qr(core,0);
    rnew=min(r3,ncur*r2);
    y{i}=reshape(y{i},[ncur,r2,rnew]);
    y{i}=permute(y{i},[1,2,3]);
   n=size(x{i},1); rx1=size(x{i},2);  rx2=size(x{i},3);
   ry1=size(y{i},2); ry2=size(y{i},3);
   %phi x is ry1 x rx1 core_x is n x rx1 x rx2
   phi{i}=phi{i-1}*reshape(permute(x{i}, [2,1,3]),[rx1,n*rx2]);
   %phi is now ry1 x n x rx2, core_y is n x ry1 x ry2  
   phi{i}=permute(reshape(phi{i},[ry1,n,rx2]),[2,1,3]); % n x ry1 x rx2
   phi{i}=reshape(y{i},[n*ry1,ry2])'*reshape(phi{i},[n*ry1,rx2]);
end
%Now we are able to recompute Y_d
%mat_x{d} is X(i_d,rx) and phi is  ry x rx
for s=1:nswp

%fprintf('Last phi-norm %e \n',norm(phi{d-1}(:)));
%fprintf('phi_true \n');
%phi{d-1}
y{d}= x{d}*phi{d-1}';
%fprintf('y_true \n');
%y{d}

%res=y;
%return;
%Now we can start right-to-left sweep that combines computation of 
%phi matrices with Y right-to-left orthogonalization
[q,rv]=qr(y{d},0);
y{d}=q;
%mat_x{d} is n x rx, mat_y{d} is n x ry phi will be ry x rx
phi{d}=y{d}'*x{d};
%keyboard;


%return
%fprintf('true phi{%d} norm: %e \n',d,norm(phi{d}(:)));

for i=d-1:-1:2
    %fprintf('i=%d/%d \n',i,d-1);
   %First solve for y{i} then orthogonalize it
   %Since the left phi matrix for y dissappears we have to
   %convolve the current core X(a1,i1,a2) with phi(i-1) from the left
   %and with phi(i) from the right
   %y{i}=ten_conv(x{i},3,phi{i+1}'); %is n x rx1 x ry2
   phi0=ten_conv(x{i},3,phi{i+1}');
   %printf('true \n'); y{i}
   %return
   
   y{i}=ten_conv(phi0,2,phi{i-1}'); %is n x ry1 x ry2
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
    
    %Compute new phi{i}
    ncur=size(x{i},1); rx1=size(x{i},2);
    ry1=size(y{i},2); ry2=size(phi{i+1},1);
    
    %X is ncur x rx1 x rx2, phi is ry2 x rx2
  % keyboard;
    %phi0=reshape(x{i},[ncur*rx1,rx2])*phi{i+1}';
    %phi0=reshape(x{i},[ncur*rx1,rx2])*phi{i+1}';
    
    %phi0 is ncur x rx1 x ry2, y{i} is ncur x ry1 x ry2
    phi0=permute(reshape(phi0,[ncur,rx1,ry2]),[1,3,2]); 
    phi{i}=reshape(permute(y{i},[1,3,2]),[ncur*ry2,ry1])'*reshape(phi0,[ncur*ry2,rx1]); 
    %fprintf('true phi{%i} norm: %e \n',i,norm(phi{i}(:)));

end
%fprintf('True: \n');
%for i=2:d
%  fprintf('%e \n',norm(phi{i}(:)));
%end
y{1} = x{1}*phi{2}';
%phi{2}
%fprintf('true phi{2} norm: %e \n',norm(phi{2}(:)));
%res=y;
%return;
%Now (finally!) we have this left-to-right sweep
[q,rv]=qr(y{1},0);
y{1}=q;
%X(i1,a)*Y(i1,b)  
phi{1}=y{1}'*x{1}; %will be ry x rx
%phi{1}
for i=2:d-1
    %fprintf('i=%d/%d \n',i,d-1);
   y{i}=ten_conv(x{i},3,phi{i+1}');
   y{i}=ten_conv(y{i},2,phi{i-1}');
   %norm(y{i}(:))
   %return
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
   n=size(x{i},1); rx1=size(x{i},2);  rx2=size(x{i},3);
   ry1=size(y{i},2); ry2=size(y{i},3);
   %phi x is ry1 x rx1 core_x is n x rx1 x rx2
   phi{i}=phi{i-1}*reshape(permute(x{i}, [2,1,3]),[rx1,n*rx2]);
   %phi is now ry1 x n x rx2, core_y is n x ry1 x ry2  
   phi{i}=permute(reshape(phi{i},[ry1,n,rx2]),[2,1,3]); % n x ry1 x rx2
   phi{i}=reshape(y{i},[n*ry1,ry2])'*reshape(phi{i},[n*ry1,rx2]);
end
%Now we are able to recompute Y_d
%mat_x{d} is X(i_d,rx) and phi is  ry x rx
end
y{d}= x{d}*phi{d-1}';
res=y;
return
end
