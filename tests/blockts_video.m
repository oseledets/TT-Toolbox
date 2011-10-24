last_indices = cell(1, d0t);

% u2 = qtt_to_tt(core(U), [d0x*ones(dpx*dconf,1); ones(d0t,1)], 0);
% u2 = tt_compr2(u2, tol);

ind = ones(1,d0t);
ind(1:5)=2;
while (max(ind<2*ones(1,d0t))==1) 
    
    for i=1:d0t
        last_indices{i}=ind(i);
    end;
%     u0 = core(U);
    u0 = tt_elem3(u2, last_indices);
    u0 = tt_squeeze(u0);
    u0 = tt_tensor(u0);
    
    n=2^d0x/2;
    contour(full(u0(n,n,n,n, :,n,n,:, n,n,n,n), 2^d0x*[1,1]));
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
    titl = sprintf('time: %g (t=%d)', t*tau, t);
%     titl = sprintf('time: %g (t=%d)', tranges(out_t)+t*tau, t);
    title([' ' titl]);
%     keyboard;
    pause(0.01);
    
    ind(6)=ind(6)+1;
    for j=1:d0t-1
        if (ind(j)>2)
            ind(j+1)=ind(j+1)+1;
            ind(j)=1;
        end;
    end;
end;