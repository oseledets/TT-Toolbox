n = 1024;
alpha = 0.01;
Nt = 1024;

a = -10;
b = 10;
T = 10;

h = (b-a)/(n-1);
tau = T/Nt;

% u = cos(pi*(a:h:b)');
u = exp(-((a:h:b).^2)');
% u = -rand(n,1);
% f = 0*ones(n,1);
f = exp(-(((a:h:b)-5).^2)') + exp(-(((a:h:b)+5).^2)');
% f = cos(pi*(a+h:h:b-h)')+1;

A = spdiags(ones(n,1)*[-1,0,1], [-1 0 1], n, n);
% A(1,1)=-1;
A(n,n)=1;
A=A/h;

M = spdiags(ones(n,1)*[1/3,1/3,1/3], [-1 0 1], n, n);

% AM = newton_prepare(A,M);

B = spdiags(ones(n,1)*[-1,2,-1], [-1 0 1], n, n);
% B(1,1)=1;
B(n,n)=1;
B = B*alpha/(h^2);

% hold on;
ttimes = zeros(Nt,1);
for t=1:Nt
    t0 = tic;
    u = burg_timestep(u, f, tau, B, A, M);
    ttimes(t)=toc(t0);
    
    max_rank = max(tt_ranks(full_to_qtt(u, 1e-10)));
    fprintf('Timestep %d (%3.5f) done. Max_rank: %d, time: %3.5f\n', t, tau*t, max_rank, ttimes(t));
    
    plot(u, '-');
    str_title = sprintf('Timestep %d (%3.5f)', t, tau*t);
    title(str_title)
    pause(0.1);
end;
% hold off;
