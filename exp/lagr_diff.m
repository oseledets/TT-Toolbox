function [D]=lagr_diff(x)

n = max(size(x));
D = zeros(n,n);

denom = ones(n,1);
for i=1:n
    for k=1:n
        if (k~=i)
            denom(i)=denom(i)*(x(i)-x(k));
        end;
    end;
end;

for j=1:n
    for i=1:n
        curD=0;
        for k1=1:n
            if (k1~=i)
                curnom = 1;
                for k2=1:n
                    if (k2~=i)&&(k2~=k1)
                        curnom=curnom*(x(j)-x(k2));
                    end;
                end;
                curD=curD+curnom;
            end;
        end;
        D(j, i)=curD/denom(i);
    end;
end;

end