function [X]=fgmres_selfprec(A, tol, restart, maxit)

n=size(A,1);
I=eye(n);

% X=I*0;
X=0*I;

% Newton:
% for i=1:30
%     RES=2*I-A*X;
%     X=X*RES;
%     err_real=norm(A*reshape(X,n,n)-I, 'fro')/norm(I, 'fro');
%     fprintf('iter [%d], res_real=%3.3e\n', i, err_real);
% end;


normf=norm(I,'fro');

H=zeros(restart+1, restart);
V=zeros(n^2, restart);
Z=zeros(n^2, restart);

for iout=1:maxit
    resid=I-A*X;
    beta=norm(resid, 'fro');
    if (iout==1)
        beta0=beta;
    end;
    V(:,1)=reshape(resid, n^2, 1)/beta;
    X_new=X;
    if (iout==1)
        X_new = I;
    end;
    
    for j=1:restart
%         if ((iout==1)&&(j==1))
%             Z(:,j)=V(:,j);
%         else
            Z(:,j)=reshape(X_new*reshape(V(:,j), n, n), n^2, 1);
%         end;
        
        w = reshape(A*reshape(Z(:,j), n, n), n^2, 1);
        for i=1:j
            H(i,j)=w'*V(:,i);
            w=w-H(i,j).*V(:,i);
        end;
        
        H(j+1,j) = sqrt((w')*w);
        if (j<restart)
            V(:,j+1)=w./H(j+1,j);
        end;
        
        [Q,R]=qr(H(1:j+1,1:j), 0);
        y=(Q(1,:)').*beta;
        y=R\y;
        
        err = norm(H(1:j+1,1:j)*y - [beta zeros(1,j)]')/beta0;
        
        X_new=X;
        for i=1:j
            X_new=X_new+y(i)*reshape(Z(:,i), n, n);
        end;
        
        err_real = norm(A*X_new-I, 'fro')/normf;
%         err_real=666;
        fprintf('iter [%d,%d], res_real=%3.3e, res=%3.3e\n', iout, j, err_real, err);
        if (err<tol)
            break;
        end;
    end;
    
%     X_new=X;
%     for i=1:j
%         X_new=X_new+y(i)*reshape(Z(:,i), n, n);
%     end;
    
    X=X_new;
    
%     mesh(X);
%     pause(2);
    
    if (err<tol)
        break;
    end;
end;

end