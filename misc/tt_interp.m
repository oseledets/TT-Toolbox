function [v]=tt_interp(u, grid, x)
% function [v]=tt_interp(u, grid, x)

d = max(size(u));
ns = tt_size(u);
v = u;

hs = cell(d,1);
for i=1:d
    hs{i} = diff(grid{i});
end;

ind = zeros(d,1);
for i=1:d
    curind=find(x(i)>grid{i}, 1, 'last');
    if (isempty(curind))
        ind(i)=1;
    else
        ind(i)=curind;
    end;
    if (ind(i)<1)
        ind(i)=1;
    end;
    if (ind(i)>=ns(i))
        ind(i)=ns(i)-1;
    end;
    v{i} = v{i}(ind(i),:,:)*(grid{i}(ind(i)+1)-x(i))/hs{i}(ind(i)) + ...
           v{i}(ind(i)+1,:,:)*(x(i)-grid{i}(ind(i)))/hs{i}(ind(i));
end;

v = tt_to_full(v);

end
