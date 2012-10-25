function [tt]=round(tt,varargin)
%Approximate QTT-Tucker with another one with specified accuracy
%   [QTT]=ROUND(QTT,EPS) Approximate QTT-Tucker tensor with relative 
%   accuracy EPS
%
%   [QTT]=ROUND(QTT,EPS,RMAX) Approximate QTT-Tucker tensor with relative
%   accuracy 
%   EPS and maximal rank RMAX. RMAX can be array of ranks or a number
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
d=tt.dphys;
core=tt.core;
tuck=tt.tuck;
eps=varargin{1};
tolcorr = 0;
for i=1:d
    tolcorr = tolcorr+tuck{i}.d;
end;
rmax = [];
if (nargin==3)
    rmax = varargin{2};
end;
ismatrix = 0;
if (isa(tuck{1}, 'tt_matrix'))
    ismatrix = 1;
    curn = cell(d,1);
    curm = cell(d,1);
end;
for i=1:d
    if (ismatrix)
        curn{i} = tuck{i}.n;
        curm{i} = tuck{i}.m;
        tuck{i} = tt_tensor(tuck{i});
    end;
   [tuck{i},rm]=qr(tuck{i},'lr');
   core{i}=ten_conv(core{i},2,rm.');
end
if (isempty(rmax))
    core=round(core,eps*sqrt(d)/sqrt(tolcorr)); 
else
    core=round(core,eps*sqrt(d)/sqrt(tolcorr),rmax); 
end;
%Round the core --- we know the result comes
%with rl orthogonality? -< No, we don't
[core, nrm] = qr(core, 'lr');
core{d} = core{d}*nrm;
rtt=rank(core); 
n=size(core);
for i=d:-1:1
   cr=reshape(core{i},[rtt(i),n(i),rtt(i+1)]);
   cr=permute(cr,[2,1,3]); cr=reshape(cr,n(i),rtt(i)*rtt(i+1));
   [u,s,v]=svd(cr,'econ');
   s=diag(s);
   r=my_chop2(s,norm(s)*eps/sqrt(tolcorr));
   u=u(:,1:r); s=s(1:r); v=v(:,1:r);
   tuck{i}=tuck{i}*(u*diag(s)); 
   if (isempty(rmax))
       tuck{i}=round(tuck{i},eps*sqrt(tuck{i}.d)/sqrt(tolcorr));
   else
       tuck{i}=round(tuck{i},eps*sqrt(tuck{i}.d)/sqrt(tolcorr), rmax);
   end;
   [tuck{i},rm]=qr(tuck{i},'lr');
   cr=rm*v';
   cr=reshape(cr,[r,rtt(i),rtt(i+1)]);
   % Shift QR to the next core block
   if (i>1)
       cr=permute(cr,[1,3,2]);
       cr = reshape(cr, r*rtt(i+1), rtt(i));
       [cr, rv] = qr(cr, 0);
       cr2 = core{i-1};
       rtuck2 = size(cr2, 2);
       cr2 = reshape(cr2, rtt(i-1)*rtuck2, rtt(i));
       cr2 = cr2*(rv.');
       rtt(i) = size(cr, 2);
       core{i-1} = reshape(cr2, rtt(i-1), rtuck2, rtt(i));
       core{i} = reshape(cr.', rtt(i), r, rtt(i+1));
   else
       core{i}=permute(cr,[2,1,3]);
   end;
end
if (ismatrix)
    for i=1:d
        tuck{i} = tt_matrix(tuck{i}, curn{i}, curm{i});
    end;
end;
tt.core=core;
tt.tuck=tuck;
return
end
