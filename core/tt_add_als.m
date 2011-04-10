function [res] = tt_add_als(x,y,nswp)
% [RES]=TT_ALS(X,Y,NSWP)
% This is a special function that computes several ALS sweeps for
% a sum of tensors in TT format. 
% X is given (cell of tensors in TT format), 
% Y is the initial guess (also in TT format)
% Output: RES --- Improved TT approximation of sum X1+...+Xm
%
% TT Toolbox 1.1, 2009-2010
%
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Contribution: Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
res=0;
m=size(x,1);
d=size(x{1},1);
%The first step is to compute left-to-right sweep of summation, yielding
%phi-matrices
phi=cell(m,1); %To store phi matrices

[q,rv]=qr(y{1},0);
y{1}=q;
for j=1:m
    phi{j}=cell(d,1);
    phi{j}{1}=y{1}'*x{j}{1};
end
       
for i=2:d-1
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
   n=size(x{1}{i},1);
   ry1=size(y{i},2);
   ry2=size(y{i},3);

   for j=1:m
       rx1=size(x{j}{i},2);  
       rx2=size(x{j}{i},3);

       phi{j}{i}=phi{j}{i-1}*...
           reshape(permute(x{j}{i},[2,1,3]),[rx1, n*rx2]);

       phi{j}{i}=permute(reshape(phi{j}{i},[ry1,n,rx2]),[2,1,3]); % n x ry1 x rx2

       phi{j}{i}=reshape(y{i},[n*ry1,ry2])'*reshape(phi{j}{i},[n*ry1,rx2]);

   end
   


end

%Now we are able to recompute Y_d
for s=1:nswp

y{d}=x{1}{d}*phi{1}{d-1}';
for j=2:m
    y{d}= y{d}+x{j}{d}*phi{j}{d-1}';
end

%Now we can start right-to-left sweep that combines computation of
%phi matrices with Y right-to-left orthogonalization
[q,rv]=qr(y{d},0);
y{d}=q;

for j=1:m
    phi{j}{d}=y{d}'*x{j}{d};
end


for i=d-1:-1:2
   %First solve for y{i} then orthogonalize it
   %Since the left phi matrix for y dissappears we have to
   %convolve the current core X(a1,i1,a2) with phi(i-1) from the left
   %and with phi(i) from the right
   phi0=cell(m,1);
   for j=1:m
       phi0{j}=ten_conv(x{j}{i},3,phi{j}{i+1}');
   end
   y{i}=ten_conv(phi0{1},2,phi{1}{i-1}');
   for j=2:m  
       y{i}=y{i}+ten_conv(phi0{j},2,phi{j}{i-1}'); %is n x ry1 x ry2
   end
   %Orthogonalize y
    ncur=size(y{i},1);
    r2=size(y{i},2);
    r3=size(y{i},3);
    core=reshape(permute(y{i},[1,3,2]),[ncur*r3,r2]);
    [y{i},rv]=qr(core,0);
    rnew=min(r2,ncur*r3);

    y{i}=permute(reshape(y{i},[ncur,r3,rnew]),[1,3,2]);
    
    %Compute new phi{i}
    ry1=size(y{i},2); 
    ry2=size(y{i},3); 
    ncur=size(x{1}{i},1); 
    tmp=reshape(permute(y{i},[1,3,2]),[ncur*ry2,ry1])';
    for j=1:m
        rx1=size(x{j}{i},2);
        ry2=size(phi{j}{i+1},1);
    
        phi0{j}=permute(reshape(phi0{j},[ncur,rx1,ry2]),[1,3,2]); 
        phi{j}{i}=tmp*reshape(phi0{j},[ncur*ry2,rx1]); 
    end

end

y{1} = x{1}{1}*phi{1}{2}';
for j=2:m
    y{1} = y{1}+x{j}{1}*phi{j}{2}';
end

%Now (finally!) we have this left-to-right sweep
[q,rv]=qr(y{1},0);
y{1}=q;
for j=1:m
    phi{j}{1}=y{1}'*x{j}{1}; 
end
%phi{1}

for i=2:d-1
   y{i}=ten_conv(x{1}{i},3,phi{1}{i+1}');
   y{i}=ten_conv(y{i},2,phi{1}{i-1}');
 
   for j=2:m
      tmp=ten_conv(x{j}{i},3,phi{j}{i+1}');
      tmp=ten_conv(tmp,2,phi{j}{i-1}');
      y{i}=y{i}+tmp;
   end
   
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
   n=size(x{1}{i},1);
   ry1=size(y{i},2); ry2=size(y{i},3);

   for j=1:m
       rx1=size(x{j}{i},2);  rx2=size(x{j}{i},3);
       phi{j}{i}=phi{j}{i-1}*reshape(permute(x{j}{i},[2,1,3]),[rx1,n*rx2]);
       phi{j}{i}=permute(reshape(phi{j}{i},[ry1,n,rx2]),[2,1,3]); % n x ry1 x rx2
       phi{j}{i}=reshape(y{i},[n*ry1,ry2])'*reshape(phi{j}{i},[n*ry1,rx2]);
   end
end
%Now we are able to recompute Y_d

end
y{d}= x{1}{d}*phi{1}{d-1}';
for j=2:m
    y{d}= y{d}+x{j}{d}*phi{j}{d-1}';
end
    
res=y;
return
end
