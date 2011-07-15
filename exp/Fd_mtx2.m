function [mat]=Fd_mtx2(a)
%[MAT]=FD_MTX2(a)
%Generates the diffusion coefficient. [0,1]^2 is supposed;

n=size(a,1)-1; %The coefficient is (n+1)x(n+1);

%The solution is (n).^2, thus the matrix will be of this size
%we need to create 4 arrays: up down left right
%
%  
%     |x y x | 
%      y O y
%     |x y x | 
%
%

%the stencil is 5-point 

ad=zeros(n,n);
  for i=1:n-1
    for j=1:n
       ad(i,j)=(a(i+1,j)+a(i+1,j+1))*0.5;
    end
  end
au=zeros(n,n);
au(2:n,1:n)=ad(1:n-1,1:n);
%for i=1:n-1
%  for j=1:n
%    au(i+1,j) = (a(i,j+1)+a(i+1,j+1))*0.5;
%  end
%end
%au(2:n,1:n-1)=ad(1:n-1,1:n-1); 
%au=ad';
al=zeros(n,n);
 for i=1:n
    for j=1:n-1
       al(i,j)=(a(i,j+1)+a(i+1,j+1))*0.5;
    end
 end
  
ar=zeros(n,n);
% for i=1:n
%    for j=1:n
%       ar(i,j)=(a(i+1,j)+a(i+1,j+1))*0.5;
%    end
%  end
%ar(1:n-1,2:n)=al(1:n-1,1:n-1); %The matrix is symmetric
ar(1:n,2:n)=al(1:n,1:n-1);
ac=zeros(n,n);
for i=1:n
  for j=1:n
    ac(i,j)=(a(i,j)+a(i,j+1)+a(i+1,j)+a(i+1,j+1));
  end
end
%2,4,3,1,5
%d
%u
%c
%l
%r
bar_a=[-al(:), -ad(:), ac(:), -au(:), -ar(:)];
mat = spdiags(bar_a, [-n, -1, 0, 1, n], n^2, n^2);
mat=mat*(n+1).^2; %Assumes [0,1]

return
end