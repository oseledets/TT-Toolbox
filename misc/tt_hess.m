function [Hess]=tt_hess(tt, indext, h)
%Numerical Hessian in a given point
%   [HESS]=TT_HESS(TT, INDEXT, H) Computes the numerical Hessian of the
%   TT-tensor TT computed at INDEXT. The grid size is supposed to be equal
%   to H
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
tt=core(tt); 
d = max(size(tt));
if (max(size(h))==1)
    h = h*ones(d,1);
end;
if (max(size(indext))==1)
    indext = indext*ones(d,1);
end;

% curind = cell(d,1);

Hess = zeros(d,d);
for i=1:d
    for j=1:d
%         for k=1:d
%             curind{k}=indext(k);
%         end;        
        curhess = tt;
        if (i==j) % Laplace_i            
%             curind{i}=indext(i)-1:indext(i)+1;
%             curhess = tt_elem3(tt, curind);
%             core=tt{i};
%             n = size(core,1); r1=size(core,2); r2=size(core,3);
%             core = reshape(core, n, r1*r2);
%             lp = spdiags(ones(n,1)*[1,-2,1], [-1,0,1], n, n)/(h(i)^2); 
%             core = lp*core;
%             core = sum(core,1)*h(i);
%             curhess{i}=reshape(core, 1, r1, r2);
            curhess{i}=(1/(h(i)^2))*(curhess{i}(indext(i)-1,:,:)-2*curhess{i}(indext(i),:,:)+curhess{i}(indext(i)+1,:,:));
            for k=1:d
                if (k~=i)
%                     curhess{k}=sum(curhess{k},1)*h(k);
                    curhess{k}=curhess{k}(indext(k),:,:);
                end;
            end;
            Hess(i,i)=full(tt_tensor(curhess));
        else % Grad_i * Grad_j
%             corei=tt{i};
%             ni = size(corei,1); r1i=size(corei,2); r2i=size(corei,3);
%             corei = reshape(corei, ni, r1i*r2i);
%             gradi = spdiags(ones(ni,1)*[-0.5,0,0.5], [-1,0,1], ni, ni)/h(i); 
%             corei = gradi*corei;
%             corei = sum(corei,1)*h(i);
%             curhess{i}=reshape(corei, 1, r1i, r2i);            
%             corej=tt{j};
%             nj = size(corej,1); r1j=size(corej,2); r2j=size(corej,3);
%             corej = reshape(corej, nj, r1j*r2j);
%             gradj = spdiags(ones(nj,1)*[-0.5,0,0.5], [-1,0,1], nj, nj)/h(j); 
%             corej = gradj*corej;
%             corej = sum(corej,1)*h(j);
%             curhess{j}=reshape(corej, 1, r1j, r2j);            
%             curind{i}=indext(i)-1:indext(i)+1;
%             curind{j}=indext(j)-1:indext(j)+1;
%             curhess = tt_elem3(tt, curind);
            curhess{i}=(0.5/h(i))*(curhess{i}(indext(i)+1,:,:)-curhess{i}(indext(i)-1,:,:));
            curhess{j}=(0.5/h(j))*(curhess{j}(indext(j)+1,:,:)-curhess{j}(indext(j)-1,:,:));            
            for k=1:d
                if (k~=i)&&(k~=j)
%                     curhess{k}=sum(curhess{k},1)*h(k);
                    curhess{k}=curhess{k}(indext(k),:,:);
                end;
            end;            
            Hess(i,j)=tt_to_full(curhess);
        end;
    end;
end;

end