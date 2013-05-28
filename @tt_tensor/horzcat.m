function [c]=horzcat(varargin)
%C=[A,B,...]
%   C=HORZCAT(A,B,...) Horizontal cat of block-row TT-tensors
%   The result has r(d+1) = sum(varargin{i}.r(d+1)). The first ranks 
%   should be equal for all arguments
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

c = varargin{1};
d = c.d;
n = c.n;
ra = c.r;
ps = c.ps;
r = ra;
for i=2:nargin
    b = varargin{i};
    if (~isempty(b))
        rb = b.r;
        if (r(1)~=rb(1))
            error('Inconsistent r(1) found, please check the arguments');
        end
        r(2:d+1) = r(2:d+1)+rb(2:d+1);
    end;
end;

cr = cell(d,1);
for i=1:d
    cr{i} = zeros(r(i), n(i), r(i+1));
    cr{i}(1:ra(i), :, 1:ra(i+1)) = reshape(c.core(ps(i):ps(i+1)-1), ra(i), n(i), ra(i+1));
end;

curr = ra;
for j=2:nargin
    b = varargin{j};
    if (~isempty(b))
        rb = b.r;
        ps = b.ps;
        cr{1}(1:r(1), :, curr(2)+1:curr(2)+rb(2)) = reshape(b.core(ps(1):ps(2)-1), rb(1), n(1), rb(2));
        for i=2:d
            cr{i}(curr(i)+1:curr(i)+rb(i), :, curr(i+1)+1:curr(i+1)+rb(i+1)) = reshape(b.core(ps(i):ps(i+1)-1), rb(i), n(i), rb(i+1));
        end;
        curr = curr+rb;
    end;
end;

c = cell2core(c, cr);

% if (isempty(b) )
%   c=a;
%   return
% elseif (isempty(a) )
%   c=b;
%   return
% end
% c=tt_tensor;
% d=a.d;
% ra=a.r;
% n=a.n;
% psa=a.ps;
% cra=a.core;
% rb=b.r;
% crb=b.core;
% psb=b.ps;
% rc=ra+rb; 
% if ( ra(1) == rb(1) )
% rc(1)=ra(1); %The only modification
% else
%   error('Incorrect sizes in [A,B], please check');    
% end
% psc=cumsum([1;n.*rc(1:d).*rc(2:d+1)]);
% crc=zeros(1,psc(d+1)-1);
% 
% %It seems: Ak
% %          Bk --- the first core 
% %all others are blocks
% 
% for i=1:d
%   cr1=cra(psa(i):psa(i+1)-1);
%   cr2=crb(psb(i):psb(i+1)-1);
%   cr1=reshape(cr1,[ra(i),n(i),ra(i+1)]);
%   cr2=reshape(cr2,[rb(i),n(i),rb(i+1)]);
%   cr3=zeros(rc(i),n(i),rc(i+1));
%   cr3(1:ra(i),:,1:ra(i+1))=cr1;
%   cr3(rc(i)-rb(i)+1:rc(i),:,rc(i+1)-rb(i+1)+1:rc(i+1))=cr2;
%   crc(psc(i):psc(i+1)-1)=cr3(:);
% end
% c.ps=psc;
% c.d=d;
% c.n=n;
% c.r=rc;
% c.core=crc;
