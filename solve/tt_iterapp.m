function y = tt_iterapp(op,afun,atype,afcnstr,x,eps,max_rank,max_swp,varargin)
%Apply TT matrix operator to TT vector and error gracefully.
%   TT_ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X) applies matrix operator AFUN to vector
%   X. If ATYPE is 'matrix, then AFUN is a matrix and the OP is applied
%   directly. OP is either 'mtimes' or 'mldivide'.
%   ATYPE and AFCNSTR are used in case of error.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X,P1,P2,...) allows extra arguments to
%   AFUN(X,P1,P2,...) although this usage is now discouraged in favor of
%   using anonymous functions.
%   AFUN(X,P1,P2,...,PN,TFLAG) should accept a TFLAG as its final input if
%   the calling function is BICG, LSQR or QMR. TFLAG is either 'transp' or
%   'notransp' depending on whether A' OP X or A OP X is required.
%   ITERAPP is designed for use by iterative methods like PCG which
%   require matrix operators AFUN representing matrices A to operate on
%   vectors X and return A*X and may also have operators MFUN representing
%   preconditioning matrices M operate on vectors X and return M\X.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2008 The MathWorks, Inc.
%   $Revision: 1.7.4.4 $ $Date: 2008/06/20 08:01:27 $


matvec_alg='mv+compr';   % produce matvec via tt_mv, and then tt_compr2
% matvec_alg='mvk3';         % produce matvec via tt_mvk2 
% matvec_alg='smv+compr'; % produce matvec via SparseMatVec
% matvec_alg='mv'; % Just a tt_mv (I guess we have to do a compression by ourselves)

if strcmp(atype,'matrix')
    switch lower(op)
        case 'mtimes'
            if (strcmp(matvec_alg,'mvk3'))
                y = tt_mvk3(afun, x, eps, 'rmax', max_rank, 'nswp', max_swp, 'verb', 0);                
            end;
            if (strcmp(matvec_alg,'mv+compr'))
                y = tt_mv(afun, x);
                y = tt_compr2(y, eps, max_rank);
            end;
            if (strcmp(matvec_alg,'smv+compr'))
                y = tt_smv(afun, x);
                y = tt_compr2(y, eps, max_rank);
            end;
            if (strcmp(matvec_alg,'mv'))
                y = tt_mv(afun, x);
            end;            
        case 'mldivide'
            if (strcmp(matvec_alg,'mvk3'))
                y = tt_mvk3(afun, x, eps, 'rmax', max_rank, 'nswp', max_swp);                
            end;
            if (strcmp(matvec_alg,'mv+compr'))
                y = tt_mv(afun, x);
                y = tt_compr2(y, eps, max_rank);
            end;
            if (strcmp(matvec_alg,'smv+compr'))
                y = tt_smv(afun, x);
                y = tt_compr2(y, eps, max_rank);
            end;
            if (strcmp(matvec_alg,'mv'))
                y = tt_mv(afun, x);
            end;                        
        otherwise
            error('MATLAB:iterapp:InvalidOp', 'Invalid operation.')
    end
else
    try
            y = afun(x,eps,max_rank,varargin{:});
    catch ME
        error('MATLAB:InvalidInput', ['user supplied %s ==> %s\n' ...
            'failed with the following error:\n\n%s'], ...
            atype,afcnstr, ME.message);
    end
end
