function [x]=als_solve_rx_2(mat, rhs, tol, maxit, x0, rx, nswp)
% function [x]=als_solve_rx_2(mat, rhs, [tol], [maxit], [x0], [rx], [nswp])

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
if (nargin<6)||(isempty(rx))
    rx=1;
end;
if (nargin<7)||(isempty(nswp))
    nswp=8;
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
    cur_x = als_solve_rx(mat, cur_rhs, tol, rx, nswp);    
%     cur_rhs = cur_rhs - tt_mat_full_vec(mat, cur_x);
    
    x = x+cur_x;
    cur_rhs = rhs - tt_mat_full_vec(mat,x);
    
    err = norm(cur_rhs)/nrmf;
    conv_fact = err_old/err;
    
    fprintf('als_solve_full: iter = %d, resid = %3.3e, rx = %d, conv_fact=%3.3f\n', i, err, rx, conv_fact);
    
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