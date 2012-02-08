function [c] = mtimes(a,b,varargin)
%C=A*B: Matrix-by-matrix, matrix-by-vector, matrix-by-number multiplication
%   [C]=MTIMES(A,B) if A is a TT-matrix and B is a TT-matrix, computes
%   matrix-by-matrix product, 
%
%   [C]=MTIMES(A,B) if A is a TT-matrix and B is TT-tensor, computes
%   matrix-by-vector product,
%
%   [C]=MTIMES(A,B) if A is a TT-matrix and B is a full tensor, computes
%   matrix-by-vector product,
%
%   [C]=MTIMES(A,B) if A is a TT-matrix and B is a number, computes their
%   product,
%
%   [C]=MTIMES(A,B) if A is a number and B is a TT-matrix, computes their
%   product
%
%
% TT Toolbox 2.1, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

if (isa(a,'tt_matrix') && isa(b,'double') && numel(b) == 1 || (isa(b,'tt_matrix') && isa(a,'double') && numel(a) == 1) )
   c=tt_matrix;
    if (  isa(b,'double') )
      c.n=a.n; c.m=a.m; c.tt=(a.tt)*b;
   else
      c.n=b.n; c.m=b.m; c.tt=a*(b.tt);
   end
elseif ( isa(a,'tt_matrix') && isa(b,'tt_tensor') && nargin == 2)
    %Put "zaglushka" in it
    %c=tt_tensor(tt_mv(core(a),core(b)));
    % return
     %Replace zaglushka
     c=tt_tensor;
     n=a.n; m=a.m; crm=a.tt.core; ttm=a.tt; psm=ttm.ps; d=ttm.d; rm=ttm.r;
     rv=b.r; crv=b.core; psv=b.ps; 
     rp=rm.*rv; 
     psp=cumsum([1;n.*rp(1:d).*rp(2:d+1)]);
     crp=zeros(psp(d+1)-1,1);
     c.d=d;
     c.r=rp;
     c.n=n;
     c.ps=psp;
     for i=1:d
        mcur=crm(psm(i):psm(i+1)-1);
        vcur=crv(psv(i):psv(i+1)-1);
        mcur=reshape(mcur,[rm(i)*n(i),m(i),rm(i+1)]);
        mcur=permute(mcur,[1,3,2]); mcur=reshape(mcur,[rm(i)*n(i)*rm(i+1),m(i)]);
        vcur=reshape(vcur,[rv(i),m(i),rv(i+1)]);
        vcur=permute(vcur,[2,1,3]); vcur=reshape(vcur,[m(i),rv(i)*rv(i+1)]);
        pcur=mcur*vcur; %pcur is now rm(i)*n(i)*rm(i+1)*rv(i)*rv(i+1)
        pcur=reshape(pcur,[rm(i),n(i),rm(i+1),rv(i),rv(i+1)]);
        pcur=permute(pcur,[1,4,2,3,5]); 
        crp(psp(i):psp(i+1)-1)=pcur(:);
     end
     c.core=crp;
     %Simple cycle through cores
    %fprintf('matrix-by-vector not implemented yet \n');
elseif ( isa(a,'tt_tensor') && isa(b,'tt_matrix') && nargin == 2)
        fprintf('vector-by-matrix not implemented yet \n');
elseif ( isa(a,'tt_matrix') && isa(b,'double') && nargin == 2 )
    %TT-matrix by full vector product
    n=a.n; m=a.m; tt=a.tt; cra=tt.core; d=tt.d; ps=tt.ps; r=tt.r;
    rb=size(b,2);
    c=reshape(b,[m',rb]);
    
    for k=1:d
      %c is rk*jk...jd*(i1..ik-1) tensor, conv over  
      %core is r(i)*n(i)*r(i+1)
      cr=cra(ps(k):ps(k+1)-1);
      cr=reshape(cr,[r(k),n(k),m(k),r(k+1)]);
      cr=permute(cr,[2,4,1,3]); cr=reshape(cr,[n(k)*r(k+1),r(k)*m(k)]);
      M=numel(c);
      c=reshape(c,[r(k)*m(k),M/(r(k)*m(k))]);
      c=cr*c; c=reshape(c,[n(k),numel(c)/n(k)]);
      c=permute(c,[2,1]);
    end
    c=c(:); c=reshape(c,[rb,numel(c)/rb]);
    c=c.';
    elseif ( isa(a,'tt_matrix') && isa(b,'tt_matrix') && nargin == 2)
    %fprintf('matrix-by-matrix not implemented yet \n');
%     c=tt_matrix(tt_mm(core(a),core(b)));
    % Time to write something working here
    c = tt_matrix;
    d = a.tt.d;
    n = a.n; m = b.m; p = a.m;
    c.n = n; c.m = m;
    tt = tt_tensor;
    ra = a.tt.r;
    rb = b.tt.r;
    tt.d = d;
    tt.r = ra.*rb;
    tt.n = n.*m;
    tt.ps = cumsum([1; tt.r(1:d).*n.*m.*tt.r(2:d+1)]);
    tt.core = zeros(tt.ps(d+1)-1, 1);
    for i=1:d
        cra = a.tt{i};
        cra = reshape(cra, ra(i), n(i), p(i), ra(i+1));
        crb = b.tt{i};
        crb = reshape(crb, rb(i), p(i), m(i), rb(i+1));
        cra = permute(cra, [1, 2, 4, 3]);
        cra = reshape(cra, ra(i)*n(i)*ra(i+1), p(i));
        crb = permute(crb, [2, 1, 3, 4]);
        crb = reshape(crb, p(i), rb(i)*m(i)*rb(i+1));
        crc = cra*crb;
        crc = reshape(crc, ra(i), n(i), ra(i+1), rb(i), m(i), rb(i+1));
        crc = permute(crc, [1, 4, 2, 5, 3, 6]);
        tt.core(tt.ps(i):tt.ps(i+1)-1) = crc(:);
    end;
    c.tt = tt;
    
elseif ( isa(a,'tt_matrix') && isa(b,'tt_tensor') && nargin > 3)
    c=mvk(a,b,varargin);
    fprintf('Krylov matrix-by-vector not implemented yet \n');
    
elseif ( isa(a,'tt_matrix') && isa(b,'tt_matrix') && nargin == 3)
    fprintf('Krylov matrix-by-matrix not implemented yet \n');

   
end
