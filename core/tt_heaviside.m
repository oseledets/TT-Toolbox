function [t]=tt_heaviside(n,L,ind)
% Heaviside vector in the TT format
% function [t]=tt_heaviside(n,L,ind)
% Returns the tt_tensor with elements 
%   theta(i-ind), i=1:N, ind=1,...,N, N=prod(n), theta(0)=1.
%
% See also: tt_unit, tt_shift
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et. al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

if (numel(n)==1)
    n = n*ones(1,L);
end;

if (numel(ind)==1)
    ind = tt_ind2sub(reshape(n, 1, []), ind);
end;

t = cell(L,1);

t{1} = zeros(1,n(1),2);
t{1}(1,:,1) = heaviside((1:n(1))-ind(1));
t{1}(1,:,2) = 1;

for i=2:L-1
    t{i} = zeros(2,n(i),2);
    t{i}(1,:,1)=[zeros(ind(i)-1,1); 1; zeros(n(i)-ind(i), 1)];
    t{i}(2,:,1)=heaviside((1:n(i))-ind(i)-1);
    t{i}(2,:,2)=1;
end;

t{L} = zeros(2,n(L),1);
t{L}(1,:,1)=[zeros(ind(L)-1,1); 1; zeros(n(L)-ind(L), 1)];
t{L}(2,:,1)=heaviside((1:n(L))-ind(L)-1);

if (L==1)
    t{1} = sum(t{1}, 1);
end;

t = cell2core(tt_tensor,t);
end

function [y]=heaviside(x)
% Heaviside function of a _real_ array
% function [y]=heaviside(x)
% Zeroes strictly negative values in x, all others are cast to 1

y = sign(x)+1;
y(y==1)=2;
y = y*0.5;
end