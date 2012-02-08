function [X]=tt_minres_selfprec(A, tol, eps, maxit, prec_type)
%Computation of the approximate TT-matrix inverse using self-prec method
%   X=TT_MINRES_SELFPREC(A,TOL,EPS,MAXIT,PREC_TYPE) Computes the
%   approximate matrix inverse of the TT-matrix A (given in TT1.0 format)
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

matvec='mm+compr';
max_rank=10;

d=max(size(A));
ns=tt_size(A);

I=tt_eye(ns, d); I=core(I); 
if (strcmp(prec_type, 'left'))
    X=tt_scal(tt_transp(A), 1e-15);
else
    X=tt_scal(I, 1e-150);
end;

normf=sqrt(tt_dot(tt_mat_to_vec(I), tt_mat_to_vec(I)));
err=1;
err_old=2;
sp=0;

for iout=1:maxit
    if (strcmp(matvec, 'mmk2'))
        resid=tt_mmk2(A, X, eps, max_rank);
    else
        resid=tt_mm(A,X);
        resid=tt_mat_compr(resid, eps,max_rank*2);
    end;
    resid=ttm_add(I, tt_scal(resid, -1));
    
    
    if (iout>1)||(strcmp(prec_type, 'left'))
        if (strcmp(matvec, 'mmk2'))
            Xresid=tt_mmk2(X, resid, eps,max_rank);
        else
            Xresid=tt_mm(X,resid);
            Xresid=tt_mat_compr(Xresid, eps,max_rank*2);
        end;
    else
        Xresid=resid;
    end;
    
    max_xrrank=max(tt_ranks(tt_mat_to_vec(Xresid)));
    
    if (strcmp(prec_type, 'left')) % we use left preconditioner (igmres)
        if (strcmp(matvec, 'mmk2'))
            w=tt_mmk2(A, Xresid, eps,max_rank);
            w=tt_mmk2(X, w, eps,max_rank);
        else
            w=tt_mm(A,Xresid);
            w=tt_mat_compr(w, eps,max_rank);
            w=tt_mm(X,w);
%             w=tt_mat_compr(w, eps);            
        end;
        
        w = tt_mat_to_vec(w);        
        beta=sqrt(tt_dot(tt_mat_to_vec(Xresid),tt_mat_to_vec(Xresid)));
        wr = tt_dot(w, tt_mat_to_vec(Xresid));
        ww = tt_dot(w,w); 
        
        H=zeros(2,1);
        H(1,1)=wr/(beta^2);
        H(2,1)=sqrt(ww/(beta^2)-wr^2/(beta^4));
        rhs=H(1,1)*beta;
        
        y = rhs/(H'*H);
        err = err*norm(H*y-[beta; 0])/beta; 
        y=y/beta;
        
%         y = tt_dot(w,w);
%         if (y~=0)
%             y = tt_dot(w, tt_mat_to_vec(Xresid))/y;
%         else
%             y=1;
%         end;                
    else % Use right preconditioner and fgmres
        if (strcmp(matvec, 'mmk2'))
            w=tt_mmk2(A, Xresid, eps,max_rank);
        else
            w=tt_mm(A,Xresid);
%             w=tt_mat_compr(w, eps);
        end;
        
        w = tt_mat_to_vec(w);        
        
        beta=sqrt(tt_dot(tt_mat_to_vec(resid),tt_mat_to_vec(resid)));
        wr = tt_dot(w, tt_mat_to_vec(resid));
        ww = tt_dot(w,w); 
        
        H=zeros(2,1);
        H(1,1)=wr/(beta^2);
        H(2,1)=sqrt(ww/(beta^2)-wr^2/(beta^4));
        rhs=H(1,1)*beta;
        
        y = rhs/(H'*H);
        err = norm(H*y-[beta; 0])/normf; 
        y=y/beta;
        
%         y2=ww;
%         if (y~=0)
%             y2 = wr/y2;
%         else
%             y2=1;
%         end;                          
    end;
    
    max_wrank=max(tt_ranks(w));
                
    Xresid=tt_mat_to_vec(Xresid);    
    
    X=tt_mat_to_vec(X);
    if (err<err_old)
        X = tt_axpy(1, X, y, Xresid, eps, max_rank);
        err_old=err;
    else
        sp=sp+1;
    end;
    max_xrank=max(tt_ranks(X));
    
    X = tt_vec_to_mat(X, ns, ns);
    
%     if (strcmp(matvec, 'mmk2'))
%         resid=tt_mmk2(A, X, eps, max_rank);
%     else
%         resid=tt_mm(A, X);
% %         resid=tt_mat_compr(resid, eps);
%     end;
%     resid=ttm_add(I, tt_scal(resid, -1));
%     resid=tt_mat_to_vec(resid);
%     
%     err_real=sqrt(tt_dot(resid,resid))/normf;    
    err_real=0;
    fprintf('iter [%d], res_real=%3.3e, x_rank=%d, Xr_rank=%d, w_rank=%d, err=%3.3e\n', iout, err_real, max_xrank, max_xrrank, max_wrank, err);    
    if (err<tol)||(sp>0)
        break;
    end;
    
end;

end