function [c]=vertcat(varargin)
%C=[A;B;...]
%   C=VERTCAT(A,B,...) Computes "addition of columns" of TT-tensors.
%   The result has r(1) = sum(varargin{i}.r(1)). The last ranks should be
%   equal for all arguments
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
    rb = b.r;
    if (r(d+1)~=rb(d+1))
        error('Inconsistent r(d+1) found, please check the arguments');
    end
    r(1:d) = r(1:d)+rb(1:d);
end;

cr = cell(d,1);
for i=1:d
    cr{i} = zeros(r(i), n(i), r(i+1));
    cr{i}(1:ra(i), :, 1:ra(i+1)) = reshape(c.core(ps(i):ps(i+1)-1), ra(i), n(i), ra(i+1));
end;

curr = ra;
for j=2:nargin
    b = varargin{j};
    rb = b.r;
    ps = b.ps;
    for i=1:d-1
        cr{i}(curr(i)+1:curr(i)+rb(i), :, curr(i+1)+1:curr(i+1)+rb(i+1)) = reshape(b.core(ps(i):ps(i+1)-1), rb(i), n(i), rb(i+1));
    end;
    cr{d}(curr(d)+1:curr(d)+rb(d), :, 1:r(d+1)) = reshape(b.core(ps(d):ps(d+1)-1), rb(d), n(d), rb(d+1));    
    curr = curr+rb;
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
% if ( ra(d+1) == rb(d+1) )
% rc(d+1)=ra(d+1); %The only modification
% else
%   error('Incorrect sizes in [A;B], please check');    
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
