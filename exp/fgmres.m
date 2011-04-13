function [sol,res,its] = fgmres (A,PRE,precfun,rhs,sol,im,maxits,tolIts) 
%-----------------------------------------------------------------------
%
% function [sol,res,its] = fgmres (A,PRE,rhs,sol) 
% restarted gmres with Krylov subspace of dim = im.  
% NOTE: this is actually fllexibe (FGMRES) -- allows
% variations in preconditioner 
%
%-----------------------------------------------------------------------
%%%%
 tolmac = 2.2e-16;
 outputG = 1;
 n = size(A,1)    ;
 its = 0    ;
% 
% main loop 
%
  while (its < maxits) 
      vv(1:n,1) = rhs - A*sol  ;
      ro = norm(vv(1:n,1),2)  ; 
      res(its+1) = ro;
      if (its  == 0) 
	tol1=tolIts*ro  ;
      end  ; 
      if (ro <= tol1 | its >= maxits)  
	return
      end
      t = 1.0/ ro   ;
      vv(1:n,1) = vv(1:n,1) * t  ;
%------------ initialize 1-st term  of rhs of hessenberg system..  ;
      rs(1) = ro  ;
%% print its/residual info 
      if (outputG) 
         fprintf(1,' its %d  res %e \n',its,ro)
      end 

      i = 0  ;
%
   while (i < im  &  (ro  >  tol1)  &  its < maxits)
      i=i+1  ;
      its = its + 1  ;
      i1 = i + 1 ; 
      z = feval(precfun,PRE,vv(:,i));
      Z(:,i) = z; 
%%                 modified GS  ;
      vv(1:n,i1) = A*z; 
      for j=1:i
           t = vv(1:n,j)'*vv(1:n,i1)  ;
           hh(j,i) = t  ;
	   vv(1:n,i1) = vv(1:n,i1) - t*vv(1:n,j)  ;
      end  ;
      t = norm(vv(1:n,i1),2)  ;
      hh(i1,i) = t  ;
      if (t  ~= 0.0)  
        t = 1.0 / t  ;
        vv(1:n,i1) = vv(1:n,i1)*t  ;
      end  ; %% IF 
%% 
    if (i ~= 1) 
%
%       perfrom previous transformations  on i-th column of h  ;
%
      for k=2:i 
         k1 = k-1  ;
         t = hh(k1,i)  ;
         hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)  ;
         hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)  ;
      end;  %% FOR 
    end  ; %% IF 
%%
      gam = sqrt(hh(i,i)^2 + hh(i1,i)^2)  ;
      if (gam  == 0.0) 
	gam = tolmac  ;
      end  ; 
%
%       determine plane rotation and update rhs of ls pb   ;
% 
      c(i) = hh(i,i)/gam  ;
      s(i) = hh(i1,i)/gam  ;
      rs(i1) = -s(i)*rs(i)  ;
      rs(i) =  c(i)*rs(i)  ;
%
%       determine res. norm. and test for convergence-  ;
%
      hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)  ;
      ro = abs(rs(i1))  ;
      sprintf(' ro = \n',ro) ; 
      res(its+1) = ro   ; 
      if (outputG) 
         fprintf(1,' its %d  res %e \n',its,ro)
      end 
   end   ;  %% end of while (im) loop 
% 
%       now compute solution. first solve upper triangular system.  ;
% 
      rs(i) = rs(i)/hh(i,i)  ;
      for  k=i-1:-1:1  ;
         t=rs(k)  ;
         for j=k+1:i  ;
            t = t-hh(k,j)*rs(j)  ;
         end  ;
         rs(k) = t/hh(k,k)  ;
      end  ;
%
%       done with back substitution..  ;
%       now form linear combination to get solution  ;
%         ;
      for j=1:i 
         sol = sol +rs(j)* Z(1:n,j)   ;
      end  ;
      if ((ro  <=  tol1) | (its >= maxits)) return; end; 
   end;  %% end while -- restart 



