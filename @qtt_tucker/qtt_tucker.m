function t = qtt_tucker(varargin)
%QTT-Tucker contructor (TT-format for the core+QTT-format for the factors)
%   T=QTT_TUCKER(TT,SZ,EPS) Splits large mode sizes of a TT-Tensor TT, if
%   SZ is a cell array with contain Q-dimensions of the factors. 
%   Example: qtt_tucker(tt, {[2,2,2], [3,3,3], [5,5]}, eps), where tt should
%   be a 2^3 x 3^3 x 5^2 tt_tensor.
%
%   T=QTT_TUCKER(QTT,SZ,EPS) Merges large several modes and computes
%   the Tucker decomposition in these long modes. In this case, sz should
%   be an ordinary row, containing the numbers of dimensions to merge
%   sz is the ordinary row, containing numbers of dimensions to merge.
%	Example: qtt_tucker(qtt, [7,7,2], eps), where qtt should be
%   7+7+2 = 16-dimensional tt_tensor, the result will be the qtt_tucker
%   with 3D core, and 7,7,2-dimensional tucker factors.
%
%   QTT_TUCKER(TT,EPS) is the same as split variant of
%   QTT_TUCKER(TT,SZ,EPS), assuming that SZ has all mode sizes equal to 2
%   (i.e., the mode sizes of the original TT-tensor are a power of 2)
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


if (nargin == 0)
    t.core=tt_tensor; %The Tucker core
    t.tuck=cell(0); %The tucker factors (in the QTT format)
    t.dphys=0; %Number of physical dimensions
    t.sz=0; %Additional dimensions for each factor
    t = class(t, 'qtt_tucker');
    return;
end

% Copy CONSTRUCTOR
if (nargin == 1) && isa(varargin{1}, 'qtt_tucker')
    t=qtt_tucker;
    t.core = varargin{1}.core;
    t.tuck = varargin{1}.tuck;
    t.sz = varargin{1}.sz;
    t.dphys=varargin{1}.dphys;
    return;
end

%from tt_tensor (ordinary one or QTT)
if ( isa(varargin{1},'tt_tensor') )
    if ( nargin == 3 )
        if (isa(varargin{2}, 'cell'))&&(isa(varargin{3}, 'double'))
            % TT with quantization
            tt=varargin{1};
            sz=varargin{2};
            eps=varargin{3};
            [fc,cr]=qtt_tucker_m(tt,sz,eps);
            t=qtt_tucker;
            t.core=cr;
            t.tuck=fc;
            t.sz=sz;
            t.dphys=numel(sz);
        else
            % QTT, varargin{2} is dims
            t = linqtt_to_qtttucker(varargin{1}, varargin{2}, varargin{3});
        end;
    else
        tt=varargin{1};
        eps=varargin{2};
        n=tt.n;
        d=tt.d;
        sz=cell(d,1);
        for i=1:d
            sz{i}=2*ones(1,log2(n(i)));
        end
        [fc,cr]=qtt_tucker_m(tt,sz,eps);
        t=qtt_tucker;
        t.core=cr;
        t.tuck=fc;
        t.sz=sz;
        t.dphys=numel(sz);        
    end;
end

% From factors and core
if (nargin==2) && isa(varargin{1}, 'cell') && isa(varargin{2}, 'tt_tensor')
    t = qtt_tucker;
    t.core = varargin{2};
    t.dphys = numel(varargin{1});
    t.tuck = cell(t.dphys, 1);
    t.sz = cell(t.dphys, 1);
    for i=1:t.dphys
        t.tuck{i} = varargin{1}{i};
        if (isa(t.tuck{i}, 'double'))
            t.tuck{i} = tt_tensor(t.tuck{i});
        end;
        t.sz{i} = t.tuck{i}.n;
    end;
%     t = class(t, 'qtt_tucker');
end;

%From full format
if (  nargin ==3 && isa(varargin{1},'double') && isa(varargin{2},'cell') && isa(varargin{3},'double'))
    %The method is as follows. %Tr
    t=qtt_tucker;
    a=varargin{1};
    sz1=numel(a);
    sz=varargin{2};
    eps=varargin{3};
    eps0=eps;
  
    %Get physical dimensions from quantized ones
    d=numel(sz);
    n=zeros(d,1);
    for i=1:d
      n(i)=prod(sz{i});
    end
    if ( sz1 ~= prod(n) )
       error('qtt_tucker: check_sizes');
    end
    rtt=zeros(d+1,1); %TT-ranks
    rtuck=zeros(d,1);
    tuck=cell(d,1); 
    rtt(1)=1; rtt(d+1)=1;
    tt_core=[];
    for k=1:d-1
       a=reshape(a,[rtt(k),n(k),numel(a)/(n(k)*rtt(k))]);
       %First, compute the 
       a=permute(a,[2,1,3]);
       a=reshape(a,[n(k),numel(a)/n(k)]);
       [u,s,v]=svd(a,'econ');
       s=diag(s);
       r=my_chop2(s,norm(s)*eps0);
       rtuck(k)=r;
       u=u(:,1:r); s=s(1:r);
       v=v(:,1:r); 
       u=u*diag(s); %This is a good  Tucker factor
       tuck{k}=tt_tensor(u,sz{k},1,r,eps0);
       [tuck{k},rm]=qr(tuck{k},'lr');
       % v is (rtt(k)*M)*r
       % v' is r*(rtt(k)*M)
       v=rm*v'; %v is in fact rxMxrtt(k)
       %v is rxrtt
       v=reshape(v,[r,rtt(k),numel(v)/(r*rtt(k))]);
       v=permute(v,[2,1,3]); 
       v=reshape(v,[rtt(k)*r,numel(v)/(r*rtt(k))]);
       [u1,s1,v1]=svd(v,'econ');
       s1=diag(s1); 
       r1=my_chop2(s1,norm(s1)*eps);
       rtt(k+1)=r1;
       u1=u1(:,1:r1);s1=s1(1:r1); v1=v1(:,1:r1);
       a=diag(s1)*v1';
       tt_core=[tt_core;u1(:)];
    end
    %Now we have to compute the last Tucker core
    a=reshape(a,[rtt(d),n(d),numel(a)/(n(d)*rtt(d))]);
    a=permute(a,[2,1,3]);
    a=reshape(a,[sz{d},numel(a)/n(d)]);
    tuck{d}=tt_tensor(a,sz{d},1,r,eps0);
    [tuck{d},rm]=qr(tuck{d},'lr');
    a=rm; a=reshape(a,[r,rtt(d),numel(a)/(r*rtt(d))]);
    a=permute(a,[1,2,3]);
    rtuck(d)=r;
    tt_core=[tt_core;a(:)];
    rtt(d+1)=size(a,3);
    cr=tt_tensor;
    cr.d=d;
    cr.n=rtuck;
    cr.r=rtt;
    cr.core=tt_core;
    cr.ps=cumsum([1;cr.n.*cr.r(1:cr.d).*cr.r(2:cr.d+1)]);
    t.core=cr;
    t.tuck=tuck;
    t.dphys=d;
    t.sz=sz;
    
end


return
