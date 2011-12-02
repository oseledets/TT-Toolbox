function [nrm] = norm(t,varargin)
%% [NRM]=NORM(TT)
% 
% Matrix norm of tt_matrix
% The behavior is similar to the Matlab built-in function.
% For tt_matrices:
% NORM(TT)/NORM(TT,'fro')/NORM(TT,'F') is the Frobenius norm of the tt_matrix
% NORM(TT,2)/NORM(TT,'2') is the 2-norm of the tt_matrix
% The 2-norm is approximated by a power iteration and is much more expensive
% than the Frobenius norm

if (nargin == 1)
  typ='F';
else
  typ=varargin{1};
end


switch typ
 case {'fro','F'} % Frobenius norm
  nrm=sqrt(abs(dot(t,t)));
  
 case {'2', 2} % 2-norm written by Thomas Mach
  eps1=1e-8;
  s=1; % oversampling p=0
  maxit = 20; % 20 iterations
  tt=t.tt;
  x = cell(s,1);
  Mx = cell(s,1);
  MU = zeros(s,s);
  for j=1:s
    r = tt_random(t.n,tt.d,1); % choose structured random test vectors
    x{j} = tt_tensor(r);
  end
% orthonormalize
  for j=1:s
    for k=1:j-1
      alpha = -dot(x{k},x{j});
      x{j} = round(x{j} + alpha*x{k} ,eps1);
    end          
    %norm(x{j})
    x{j} = round(x{j}/norm(x{j}) ,eps1);
  end;
  for i=1:maxit
    for j=1:s
%      Mx{j}=round(t*x{j},eps1);
      x{j}=round(t*x{j},eps1);
      x{j}=round(t'*x{j},eps1);
    end

% % diagonalizing
%     for j=1:s
%       for k=1:s
% 	MU(j,k)=dot(x{k},Mx{j});
%       end
%     end
%     [V,D]=eig(MU);
%     [~,I]=sort(abs(diag(D)),'descend');
%     V=V(:,I);
%     %diag(D(I,I))
%     
%     for j=1:s
%       x{j}=0*x{j};
%     end
%     for j=1:s
%       for k=1:s
% 	x{k} = round(x{k}+Mx{j}*V(j,k),eps1);
%       end            
%     end

    if (s>1)
      s=s-1;
      MU=MU(1:s,1:s);
    end
 
    % orthonormalize
    for j=1:s
      for k=1:j-1
	alpha = -dot(x{k},x{j});
	x{j} = round(x{j} + alpha*x{k} ,eps1);
      end          
      x{j} = round(x{j}/norm(x{j}) ,eps1);
    end; 
  end % Iteration

  for j=1:s
    Mx{j}=round(t*x{j},eps1);
    Mx{j}=round(t'*Mx{j},eps1);
    for k=1:s
      MU(j,k)=dot(x{k},Mx{j});
    end
  end
  nrm=sqrt(max(eig(MU))); % approx. from below to the smallest eigenvalue
return
end