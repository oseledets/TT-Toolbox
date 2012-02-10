function [x]=als_solve_rx_2(mat, rhs, tol, maxit, x0, drx, nswp)
%Computes an approximate low-rank solution for 2D case (Method 2)
%   [x]=ALS_SOLVE_RX(MAT, RHS, [TOL], [MAXIT],[X0],[DRX], [NSWP])
%   Finds a solution to 2D TTM matrix MAT using the ALS to a 2D TT tensor
%   with rank rx, but the RHS and X are represented as full vectors.
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


nrmf = norm(rhs);

if (nargin<3)||(isempty(tol))
    tol=1e-12;
end;
if (nargin<4)||(isempty(maxit))
    maxit=2;
end;
if (nargin<5)||(isempty(x0))
    x0 = zeros(size(rhs));
end;
if (nargin<6)||(isempty(drx))
    drx=1;
end;
if (nargin<7)||(isempty(nswp))
    nswp=4;
end;


x=x0;
cur_rhs = rhs - tt_mat_full_vec(mat,x);

err = norm(cur_rhs)/nrmf;
if (err<tol) 
    return;
end;

spunct = 0;
err_old = err;
for i=1:maxit
    cur_x = als_solve_rx(mat, cur_rhs, tol, drx, nswp);    
%     cur_rhs = cur_rhs - tt_mat_full_vec(mat, cur_x);
    
    x = x+cur_x;
    cur_rhs = rhs - tt_mat_full_vec(mat,x);
    
    err = norm(cur_rhs)/nrmf;
    conv_fact = err_old/err;
    
    fprintf('als_solve_full: iter = %d, resid = %3.3e, conv_fact=%3.3f\n', i, err, conv_fact);
    
    if (conv_fact<1.2)
        spunct = spunct+1;
%         if (rx<15)
%             rx=rx+1;
%         end;
    end;
    if (spunct>3)
        break; % shit happened - we've stagnated
    end;
    err_old = err;
    if (err<tol)
        break;
    end;

%     if (rx==10)
%         break;
%     end;
%     if (conv_fact<1.01)
%         break;
%     end;
end;

end
