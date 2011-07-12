function [b]=convten(a,x,mode)
% function [b]=convten(a,x,mode)
% Convolves the tt_array A with the tt_tensor X over the mode MODE
% The order of modes in the result is the same except the convolved one

d = x.d;
rx = x.r;
psx = x.ps;
crx = x.core;
nx = x.n;

nsa = a.ns;
dp = size(nsa,2);
tta = a.tt;
ra = tta.r;
psa = tta.ps;
cra = tta.core;

% We need to determine, which modes left
nsb = zeros(d, dp-1);
cnt = 1;
% For the permutations later
perm = zeros(1,dp);
for i=1:dp
    if (i~=mode)
        nsb(:,cnt)=nsa(:,i);
        perm(cnt)=i;        
        cnt=cnt+1;
    end;
end;
perm(dp)=mode;
perm = perm+1;

rb = rx.*ra;
psb = cumsum([1; rb(1:d).*prod(nsb,2).*rb(2:d+1)]);
% Preallocate
crb = zeros(psb(d+1),1);

for i=1:d
    a1 = cra(psa(i):psa(i+1)-1);
    a1 = reshape(a1, [ra(i), nsa(i,:), ra(i+1)]);
    x1 = crx(psx(i):psx(i+1)-1);
    x1 = reshape(x1, [rx(i), nx(i), rx(i+1)]);
    
    x1 = permute(x1, [2 1 3]);
    x1 = reshape(x1, nx(i), rx(i)*rx(i+1));
    % X now is ready
    
    % Permute so that the convolving mode is the last one
    a1 = permute(a1, [1,dp+2,perm]);    
    a1 = reshape(a1, ra(i)*ra(i+1)*prod(nsb(i,:),2), nsa(i,mode));
    
    % Multiply
    b1 = a1*x1; % size ra1*ra2*nsb*rx1*rx2
    % Permute in the form [ra*rx, nsb, ra*rx]
    b1 = reshape(b1, [ra(i),ra(i+1),nsb(i,:),rx(i),rx(i+1)]);
    b1 = permute(b1, [1,3+dp-1,(3:3+dp-2),2,4+dp-1]);
    crb(psb(i):psb(i+1)-1) = reshape(b1, rb(i)*prod(nsb(i,:),2)*rb(i+1), 1);
end;

ttb = tt_tensor;
ttb.core = crb;
ttb.d = d;
ttb.n = prod(nsb,2);
ttb.r = rb;
ttb.ps = psb;

b = tt_array;
b.ns = nsb;
b.tt = ttb;

end