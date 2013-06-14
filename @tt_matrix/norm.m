function [nrm] = norm(t,varargin)
%Matrix norm of the TT-matrix
%   The behavior is similar to the Matlab built-in function.
%   NRM=NORM(TT)/NORM(TT,'fro')/NORM(TT,'F') is the Frobenius norm of the
%   TT-matrix
%   
%   NRM=NORM(TT,2)/NORM(TT,'2') is the 2-norm of the TT-matrix. The 2-norm 
%   is approximated by a power iteration and is much more expensive than 
%   the Frobenius norm
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
if (nargin == 1)
  typ='F';
else
  typ=varargin{1};
end


switch typ
 case {'fro','F'} % Frobenius norm
  %nrm=sqrt(abs(dot(t.tt,t.tt, true)));
  nrm = norm(t.tt); %More robust
 case {'2', 2} % 2-norm written by Thomas Mach
  maxit = 20; % 20 iterations
  tol = 0.01; % relative norm change
  eps1=1e-8; % tt-approximation accuracy
  rmax = 100; % maximal rank
  diag = 0; % no diagonalization
  s=1; % oversampling p=0
  for i=2:2:length(varargin)-1
    switch lower(varargin{i})
     case 'maxit'
      maxit=varargin{i+1};
     case 'tol'
      tol=lower(varargin{i+1});
     case 'x0'
      x0=varargin{i+1};
     case 'eps'
      eps1=varargin{i+1};
     case 'rmax'
      rmax=varargin{i+1};
     case 'diag'
      diag=varargin{i+1};
     case 's'
      s=varargin{i+1};
     otherwise
      error('Unrecognized option: %s\n',varargin{i});
    end
  end  
  tt=t.tt;
  x = cell(s,1);
  Mx = cell(s,1);
  MU = zeros(s,s);
  
  if (exist('x0','var'))
    x{1}=x0;
  else
    r = tt_random(t.n,tt.d,1); % choose structured random test vectors
    x{1} = tt_tensor(r);
  end    
  for j=2:s
    r = tt_random(t.n,tt.d,1); % choose structured random test vectors
    x{j} = tt_tensor(r);
  end
% orthonormalize
  for j=1:s
    for k=1:j-1
      alpha = -dot(x{k},x{j});
      x{j} = round(x{j} + alpha*x{k} ,eps1,rmax);
    end          
    %norm(x{j})
    x{j} = round(x{j}/norm(x{j}) ,eps1,rmax);
  end;
  n2 = inf;
  for i=1:maxit
    for j=1:s
      if (diag)
	Mx{j} = round(t*x{j},eps1,rmax);
      end
      x{j} = round(t*x{j},eps1,rmax);
      x{j} = round(t'*x{j},eps1,rmax);
    end

% % diagonalizing
    if (diag)
      for j=1:s
	for k=1:s
	  MU(j,k)=dot(x{k},Mx{j});
	end
      end
      [V,D]=eig(MU);
      [~,I]=sort(abs(diag(D)),'descend');
      V=V(:,I);
      %diag(D(I,I))
      
      for j=1:s
	x{j}=0*x{j};
      end
      for j=1:s
	for k=1:s
	  x{k} = round(x{k}+Mx{j}*V(j,k),eps1,rmax);
	end            
      end
    end
    
% orthonormalize
    n1 = norm(x{1});
    x{1} = round(x{1}/n1,eps1,rmax);
    for j=2:s
      for k=1:j-1
	alpha = -dot(x{k},x{j});
	x{j} = round(x{j} + alpha*x{k} ,eps1,rmax);
      end          
      x{j} = round(x{j}/norm(x{j}) ,eps1,rmax);
    end
    if (abs(n2-n1)/abs(n1)<tol) break; % termination by tolerance
    end 
    n2 = n1;
  end % Iteration
  
  for j=1:s
    Mx{j}=round(t * x{j},eps1,rmax);
    Mx{j}=round(t' * Mx{j},eps1,rmax);
    for k=1:s
      MU(j,k) = dot(x{k},Mx{j});
    end
  end
  nrm = abs(sqrt(max(eig(MU)))); % approx. from below to the largest singular value
  %The norm should be real 
  return
 otherwise
     warning('unknown norm type');
     nrm=-1;
end
end
