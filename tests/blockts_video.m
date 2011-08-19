last_indices = cell(1, d0t);

ind = ones(1,d0t);
ind(1:0)=2;
while (max(ind<2*ones(1,d0t))==1) 
    
    for i=1:d0t
        last_indices{i}=ind(i);
    end;
    u0 = core(U);
    u0 = tt_elem3(u0, last_indices);
    u0 = tt_squeeze(u0);
    u0 = tt_tensor(u0);
    
    if (dphys==1)
        plot(x,full(u0,2^d0x));
    else
%         mesh(full(u0, [2^d0x, 2^d0x]));
        scatter3(x,y,full(u0,4^d0x),5,full(u0,4^d0x));
        colorbar;
        view(2);
    end;
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
    
    ind(1)=ind(1)+1;
    for j=1:d0t-1
        if (ind(j)>2)
            ind(j+1)=ind(j+1)+1;
            ind(j)=1;
        end;
    end;
end;