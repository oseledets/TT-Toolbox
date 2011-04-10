function [qtt]=qtt_newton(gridx, gridy, d, param)
% [qtt]=qtt_newton(gridx, gridy, d, param)
% param in {1...10 } the larger, the better accuracy

eps=1e-10;
n=size(gridx,1)-1;

switch param
   case 1
       a=-15;   b=10;   rmax=65;
   case 2
       a=-15;   b=10;   rmax=80;
   case 3
       a=-20;   b=15;   rmax=125;
   case 4
       a=-25;   b=20;   rmax=145;
   case 5
       a=-25;   b=20;   rmax=200;
   case 6
       a=-25;   b=20;   rmax=220;
   case 7
       a=-25;   b=20;   rmax=245;
   case 8
       a=-30;   b=25;   rmax=320;
   case 9
       a=-35;   b=30;   rmax=410;
   case 10
       a=-35;   b=30;   rmax=440;
    otherwise
       a=-35;   b=30;   rmax=475;
end

h=(b-a)/(rmax-1);
s=(1:rmax);     s=(s-1)*h+a;

res=cell(3,1);
for i=1:3
    res{i}=zeros(n, n, rmax);
end

w=zeros(rmax,1);
for alpha=1:rmax
    w(alpha)=h*exp(s(alpha))*2/sqrt(pi);
end
w(1)=0.5*w(1);   w(rmax)=0.5*w(rmax);

for alpha=1:rmax
    for dim=1:3 
        for j=1:n
            res{dim}(j,:,alpha) = ...
                (erf((gridx(j+1,dim)-gridy(:,dim))*exp(s(alpha)))-...
                erf((gridx(j,dim)-gridy(:,dim))*exp(s(alpha))))*...
                sqrt(pi)*0.5/exp(s(alpha));
        end
    end
end
for alpha=1:rmax
    res{1}(:,:,alpha)=res{1}(:,:,alpha)*w(alpha);
end

tmp=cell(3,1);
for dim=1:3
   tmp{dim}=reshape(res{dim}(:,:,1), 2*ones(1, 2*d));
   tmp{dim}=tt_matrix(tmp{dim}, eps);
end
summand=kron(tmp{1},tmp{2});
summand=kron(summand,tmp{3});
qtt=summand;

for alpha=2:rmax
    for dim=1:3
        tmp{dim}=reshape(res{dim}(:,:,alpha), 2*ones(1, 2*d));
        tmp{dim}=tt_matrix(tmp{dim}, eps);
    end
    summand=kron(tmp{1},tmp{2});
    summand=kron(summand,tmp{3});
    qtt=qtt+summand;
    qtt=round(qtt, eps);
end

end


