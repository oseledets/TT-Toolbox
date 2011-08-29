function [p]=lagr_int(xb, xc)

n = max(size(xb));
M = max(size(xc));
p = zeros(M, n);

for j=1:M
    for i=1:n
        curp=1;
        for k=1:n
            if (k~=i)
                curp=curp*(xc(j)-xb(k))/(xb(i)-xb(k));
            end;
        end;
        p(j, i)=curp;
    end;
end;

end