rng = 1;
tdig = 4;
viewdims = [3,4];

d0t = d0ts(rng);
tau = (Trange(rng+1)-Trange(rng))/(2^d0t);

last_indices = cell(1, d0t);

% u2 = qtt_to_tt(core(Us{rng}), [d0x*ones(dpx*dconf,1); ones(d0t,1)], 0);
uft = tt_reshape(Us{rng}, [2^d0x*ones(1,dpx*dconf), 2*ones(1,d0t)]);
% u2 = tt_compr2(u2, tol);

ind = ones(1,d0t);
ind(1:tdig-1)=2;
while (max(ind<2*ones(1,d0t))==1) 
    
    cind = num2cell([1:2^d0x]'*ones(1,dpx*dconf), 1);
    cind(dpx*dconf+1:dpx*dconf+d0t) = num2cell(ind, 1);
    uf = uft(cind);
    uf = tt_reshape(uf, 2^d0x*ones(dpx*dconf, 1));
    
    n=2^d0x/2+1;
%     plin = (2^(d0x-1)-2^(d0x-2)):2^(d0x-8):(2^(d0x-1)+2^(d0x-2));
    plin=1:8:2^d0x;
    cind = num2cell(n*ones(1, dpx*dconf), 1);
    cind(viewdims(1))={plin};
    cind(viewdims(2))={plin};
%     plin = 1:max(2^(d0x-7),1):2^d0x;
    contour(full(uf(cind), max(size(plin))*[1,1]));
    view(2);
%     if (dphys==1)
%         plot(x,full(u0,2^d0x));
%     else
% %         mesh(full(u0, [2^d0x, 2^d0x]));
%         scatter3(x,y,full(u0,4^d0x),5,full(u0,4^d0x));
%         colorbar;
%         view(2);
%     end;
    t = 0;
    for j=1:d0t
        t=t+(ind(j)-1)*(2^(j-1));
    end;
    t=t+1;
    titl = sprintf('time: %g (t=%d), dimensions: [%d,%d]', Trange(rng)+t*tau, t, viewdims(1), viewdims(2));
%     titl = sprintf('time: %g (t=%d)', tranges(out_t)+t*tau, t);
    title([' ' titl]);
%     keyboard;
    pause(0.01);
    
    ind(tdig)=ind(tdig)+1;
    for j=1:d0t-1
        if (ind(j)>2)
            ind(j+1)=ind(j+1)+1;
            ind(j)=1;
        end;
    end;
end;