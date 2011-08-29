function [coefs,nodes,cw]=weight_polys(n, wfun)

coefs = zeros(n+1,n+1);
[nodes, cw]=fclencurt(n+1, -1, 1);
% nodes = cos((2*(1:n+1)'-1)*pi*0.5/(n+1));

for i=0:n % poly degree    
    coefs(i+1,i+1)=1; % the current degree
    for j=0:i-1 % orthogonalize to the previous ones
        dotprod = 0;
        for k=1:n+1 % Clenshau-Kurtis quadrature
            Pj = 0;
            for m=j:-1:0 % current poly
%                 Pj = Pj + coefs(j+1,m+1)*(nodes(k)^m);
                Pj = Pj + coefs(j+1,m+1)*cos(m*acos(nodes(k)));
            end;
            Pi = 0;
            for m=i:-1:0 % current poly
%                 Pi = Pi + coefs(i+1,m+1)*(nodes(k)^m);
                Pi = Pi + coefs(i+1,m+1)*cos(m*acos(nodes(k)));
            end;
            dotprod = dotprod + Pi*Pj*wfun(nodes(k))*cw(k);
        end;
        % remove the projection on the previous poly
        for m=0:j
            coefs(i+1,m+1)=coefs(i+1,m+1)-coefs(j+1,m+1)*dotprod;
        end;
    end;
    % Normalize
    dotprod = 0;
    for k=1:n+1 % Clenshau-Kurtis quadrature
        Pi = 0;
        for m=i:-1:0 % current poly
%             Pi = Pi + coefs(i+1,m+1)*(nodes(k)^m);
            Pi = Pi + coefs(i+1,m+1)*cos(m*acos(nodes(k)));
        end;
        dotprod = dotprod + (Pi^2)*wfun(nodes(k))*cw(k);
    end;
    for m=i:-1:0
        coefs(i+1,m+1)=coefs(i+1,m+1)/sqrt(dotprod);
    end;
end;

end