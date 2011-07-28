% Uses Ag, Pg, fg - matrix, Prec and rhs, solution is ug

maxit = 50;
eps = 1e-6;
tol = eps;
maxrank=[];

% ug=fg;
% ug = tt_tensor(tt_random(fg.n, fg.d, 2));

if (~isempty(Pg))
    Pfg = mvk2(Pg, fg, eps, 10);
else
    Pfg = fg;
end;
normfg = norm(Pfg);

ds_results = zeros(maxit,5);
for i=1:maxit
    tstart = tic;
%     ug = dmrg_solve2(Ag,fg,ug,eps,tol,maxrank,1,Pg,false);
    ug = dmrg_solve2(Ag,fg,eps,'x0',ug,'nswp',1,'P',Pg,'verb',false);
    cur_time = toc(tstart);
    
    resg = mvk2(Ag,ug,eps,10);
    if (~isempty(Pg))
        resg = mvk2(Pg, resg, eps, 10);
    end;
    resid = norm(resg-Pfg)/normfg;
    
    ds_results(i,1)=i;
    ds_results(i,2)=resid;
    ds_results(i,4)=cur_time;
    ds_results(i,3)=erank(ug);
    
    ds_results(i,5)=sum(ds_results(1:i,4));
    
    fprintf('\nsweep\t resid \t\t erank  \t sw. time  \t full time\n');
    fprintf('%d\t %3.3e\t %3.3f  \t %3.4f  \t %3.4f\n', ds_results(i,1),ds_results(i,2),ds_results(i,3),ds_results(i,4),ds_results(i,5));

    if (resid<tol) 
        break;
    end;
end;


fprintf('sweep\t resid \t erank  \t sw. time \t full time\n');
for i=1:maxit
    fprintf('%d\t %3.3e\t %3.3f  \t %3.4f \t %3.4f\n', ds_results(i,1),ds_results(i,2),ds_results(i,3),ds_results(i,4),ds_results(i,5));
end;
