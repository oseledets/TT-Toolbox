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
    
    mesh(full(u0, [2^d0x,2^d0x]));
    t = 0;
    for j=1:d0t
        t=t+(ind(j)-1)*(2^(j-1));
    end;
    t=t+1;
    titl = sprintf('time: %g (t=%d)', tranges(out_t)+t*tau, t);
    title([' ' titl]);
%     keyboard;
    pause(0.1);
    
    ind(1)=ind(1)+1;
    for j=1:d0t-1
        if (ind(j)>2)
            ind(j+1)=ind(j+1)+1;
            ind(j)=1;
        end;
    end;
end;