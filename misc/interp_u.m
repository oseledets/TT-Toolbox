function [v]=interp_u(ind,u,X,gridold,dpx,dconf)
% d = max(size(u));
d0x = log2(size(u{1},1));
d = d0x*dpx*dconf;
% d0x = round(d/(dpx*dconf));

curx = zeros(dpx*dconf,1);
indcell = cell(d,1);
for i=1:d
    indcell{i}=ind(i);
end;

for j=1:dpx
    for i=1:dconf
%         curxnew = tt_to_full(tt_elem3(X{i,j},indcell));
        curxnew = full(X{i,j}(indcell));
        if (max(size(curxnew))>1)
            keyboard;
        end;
        curx(i+(j-1)*dconf)=curxnew;
    end;
end;

% u = qtt_to_tt(u, d0x*ones(1,dpx*dconf), 0);
v = tt_interp(u,gridold,curx);
end