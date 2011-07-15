last_indices = cell(1, dt);

ind = ones(1,dt);
ind(1:5)=2;
while (sum(ind)<=2*dt) 
    
    for i=1:dt
        last_indices{i}=ind(i);
    end;
    u0 = core(U);
    u0 = tt_elem3(u0, last_indices);
    u0 = tt_squeeze(u0);
    u0 = tt_tensor(u0);
    
    plot(full(u0, 2^dx));
%     mesh(full(u0, [2^d0x,2^d0x]));
    t = 0;
    for j=1:dt
        t=t+(ind(j)-1)*(2^(j-1));
    end;
    t=t+1;
    titl = sprintf('time: %g (t=%d)', t*tau, t);
    title([' ' titl]);
%     keyboard;
    pause(0.1);
    
    ind(6)=ind(6)+1;
    for j=1:dt-1
        if (ind(j)>2)
            ind(j+1)=ind(j+1)+1;
            ind(j)=1;
        end;
    end;
end;