function [x] = amen_block_solve(A, y, tol, varargin)
%   Solution of a block linear system in the TT format by the AMEn iteration
% function [x] = amen_block_solve(A, y, tol, varargin)
% Tries to solve Ax=y, where y = [y{1}; ...; y{K}], and
% A = [A{1,1} ... A{1,K}
%       ...   ...  ...
%      A{K,1} ... A{K,K}];
% Cells A{i,j} contain tt_matrix or TT-canonical ("{d,R}") formats.
% Some cells of A can be empty, which means that the corresponding
% submatrices are zeros. All cells of y must be tt_tensors of the same
% dimensions.
%
% x is kept in the block TT format, x(i) = x^1(i_1)...x^q(i_q,k)...x^d(i_d)
% That is, the block component enumerator k=1,...,K moves through the
% train. However, the components of A and y are kept as separate TT formats.
% The solution is returned with k being either the first or the last rank;
% see the parameter 'exitdir' below.
%
% Optional parameters are given in varargin in the form 
% 'param_name1', param_value1, 'param_name2', param_value2, and so on.
% Available parameters and their default values are as follows.
%   o nswp:     number of AMEn sweeps [20].
%   o x0:       initial guess [random rank-K tensor]
%   o kickrank: residual and enrichment rank [4]
%   o max_full_size: maximal size of local systems which is treated by "\" [0]
%   o local_iters:  maximal number of local gmres iterations [1000]
%   o resid_damp:   if the local system is solved iteratively, its accuracy
%                   is tol/sqrt(d-1)/resid_damp [100]
%   o tol_exit: the method stops if all relative local residuals in one 
%               sweep are less than tol_exit [tol]
%   o exitdir:  if 1: return the solution with k in the last block, x^d(i_d,k),
%               if -1: return the solution with x^1(k,i_1) [1]
%   o auxmat:   auxiliary matrices (see below), which are to be projected 
%               onto the ALS bases []
%   o solfun:   User-defined function for solution of local systems []
% It should be of the form fun(q, XAX1,Ai,XAX2, XY1,yi,XY2, tol, sol_prev), where
% q is the TT-block number, XAX1,Ai,XAX2 are KxK cell arrays defining the
% reduced TT formats of the submatrices of the local matrix.
% XY1,yi,XY2 form the RHS vector (in the same structure as y), 
% tol is the solution tolerance, sol_prev is the vector of the previous solution. 
% If solfun returns void, the default GMRES or backslash method will be used.
% For q==1, the first TT rank is 1, and XAX1,XY1 always contain ones.
% For q==d, the last TT rank is 1, so XAX2,XY2 contain ones.
%
% If the initial matrix was a sum of kron products of sparse blocks (a 2D
% cell array), then
% numel(Ai{i,j})==numel(XAX1{i,j})==numel(XAX2{i,j}) is the number of
% Kronecker summands in A{i,j}.
% For each k, size(XAX1{i,j}{k})==[rxi1,rxj1,1],
%             size(Ai{i,j}{k}) == [ni,nj] (since Ai{i,j}{k} is sparse),
%             size(XAX2{i,j}{k})==[1,rxi2,rxj2].
%
% If A{i,j} was given as a tt_matrix, then all numels above are ones, and
%             size(XAX1{i,j}{1})==[rxi1,rxj1,ra1],
%             size(Ai{i,j}{1}) == [ra1,ni,nj,ra2],
%             size(XAX2{i,j}{1})==[ra2,rxi2,rxj2].
%
% The right-hand side blocks satisfy
%             size(XY1{i}) == [rx1,ry1],
%             size(yi{i})  == [ry1,n,ry2],
%             size(XY2{i}) == [ry2,rx2].
%
% Normally, solfun should return either the solution vector of size rx1*n*rx2*K,
% or void, indicating that we haven't computed this block and delegate this
% task to default solver. However, the return can also be cell(1,2) in the
% following situation.
% At one of TT blocks, the system can be expanded, e.g.
% A = [A{1,1}...A{1,K}  B{1}
%       ...      ...    ...
%      A{K,1}...A{K,K}  B{K}
%       ...      ...    ...
%      C{1}  ... C{K}    D  ];
% This is useful for Stokes-type equations: the matrix for the "spatial"
% TT block must be [A,B'; B, zeros], but
% on "time-stochastic" blocks we can use only the velocity part A.
% In other words, the second component in X is updated only in the spatial
% TT block, but kept freezed in the remaining TT blocks.
% However, to make the right-hand side consistent, we have to subtract the
% terms B{k}*x{K+1} from y{k} when we leave the spatial block.
% The submatrices B{k} are user-defined. However, we have to compute
% ALS reductions of their other-than-spatial TT blocks. To perform this,
% the auxiliary matrices B{1:K} must be given in the parameter 'auxmat'
% in the form cell(K,1), where auxmat{k} contains matrix TT format, either
% in tt_matrix or {d,R} form.
% In this case, solfun is obligatory and should have the form
% fun(q, XAX1,Ai,XAX2, XY1,yi,XY2, tol, sol_prev, XBX1,Bi,XBX2), where the
% last three inputs carry the necessary reductions of auxmat.
% XBX1,Bi,XBX2 are all cell(K,1), with e.g. XBX1{k} is a cell(1,R), and
% XBX1{k}{m} is a three-dimensinal array, in the same form as for A.
% However, one of TT blocks of B can be rectangular. It is the block, at
% which the system is expanded. When solfun solves this expanded system,
% its result should be cell(1,2), where result{1} contains the _active_
% (first K) solution components in the usual form rx1*n*rx2*K, while
% result{2} of size rx1*m*rx2 is the extra component x{K+1} which is frozen until
% the next iteration.

nswp = 20;
x = [];
kickrank = 4;
local_iters = 1000; % only for gmres
max_full_size = 0; % don't use the backslash by default
resid_damp = 1e2;
% Exit tolerance
tol_exit = tol;
exitdir = 1; % Sweep direction in which the method should exit
solfun = []; % Solution function
auxmat = []; % Auxiliary matrices.

% Parse parameters
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'x0'
            x=varargin{i+1};
        case 'blocking'
            warning('blocking parameter is deprecated. X is stored in a single block TT');
        case 'kickrank'
            kickrank = varargin{i+1};
        case 'solfun'
            solfun = varargin{i+1};
        case 'local_iters'
            local_iters = varargin{i+1};
        case 'max_full_size'
            max_full_size = varargin{i+1};            
        case 'resid_damp'
            resid_damp = varargin{i+1};
        case 'tol_exit'
            tol_exit = varargin{i+1};            
        case 'exitdir'
            exitdir = varargin{i+1};
        case 'auxmat'
            auxmat = varargin{i+1};
            
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end;
end;

if (isempty(solfun))
    if (~isempty(auxmat))
        error('With auxiliary matrices, the user should specify a custom solfun');
    end;
    % Default solution function with GMRES, if user did not give something else
    solfun = @(i, XAX1,Ai,XAX2, XY1,yi,XY2, tol, sol_prev)loc_solve_default(XAX1, Ai, XAX2, XY1,yi,XY2, tol, sol_prev, local_iters, max_full_size);
end;

% Block size
K = numel(y);
if (size(A,1)~=K)||(size(A,2)~=K)
    error('Block sizes of A and y differ');
end;

% Extract sizes and TT blocks of A and y
d = [];
ra = cell(K,K);
ry = cell(K,1);
for kj=1:K
    if (~isempty(y{kj}))
        if (isempty(d))
            % This is the first RHS part, determine the sizes
            d = y{kj}.d;
            n = y{kj}.n;
        else
            % All sizes must be equal
            if (y{kj}.d~=d)||(any(y{kj}.n~=n))
                error('Mismatching TT dimensionalities in RHS (%d)', kj);
            end;
        end;
        % Copy TT ranks and blocks
        ry{kj} = y{kj}.r;
        y{kj} = core2cell(y{kj});
        % y is stored as double cell y{kj}{i}, kj=1,...,K, i=1,...,d.
    else
        error('Empty right hand side blocks are not yet supported');
    end;
    
    for ki=1:K
        if (~isempty(A{ki,kj}))
            if (isa(A{ki,kj}, 'tt_matrix'))
                % Copy TT ranks and blocks from a tt_matrix
                ra{ki,kj} = A{ki,kj}.r;
                if (A{ki,kj}.d~=d)||(any(A{ki,kj}.n~=n))||(any(A{ki,kj}.m~=n))
                    error('Mismatching TT dimensionalities in matrix (%d,%d)', ki, kj);
                end;
                A{ki,kj} = core2cell(A{ki,kj});
                % A is always stored as double cell, A{ki,kj}{i,j}, where
                % ki,kj=1,...,K, i=1,...,d and j=1,...,R, where R is the
                % "canonical" rank. R==1 if a single tt_matrix is given.
                % The rank of A is a "cell-numeric" array, ra{ki,kj}(i,j).
            else
                % The matrix is given as a {d,R} cell array, meaning a sum
                % of R TT formats (needed for sparse blocks, in which case
                % it's the sum of R Kronecker products)
                R = size(A{ki,kj},2);
                if (size(A{ki,kj},1)~=d)
                    error('Mismatching TT dimensionalities in matrix (%d,%d)', ki, kj);
                end;
                ra{ki,kj}=ones(d+1,R);
                for j=1:R
                    for i=1:d
                        if (issparse(A{ki,kj}{i,j}))
                            % Sparse matrix can only have TT ranks 1
                            [n1,n2]=size(A{ki,kj}{i,j});
                            ri = 1; ri2 = 1;
                        else
                            [ri,n1,n2,ri2]=size(A{ki,kj}{i,j});
                        end;
                        if (n1~=n(i))||(n2~=n(i))
                            error('Mismatching TT dimensionalities in matrix (%d,%d)', ki, kj);
                        end;
                        ra{ki,kj}(i,j) = ri;
                        ra{ki,kj}(i+1,j) = ri2;
                    end;
                end;
            end;
        end;
    end;
end;

% Create X.
if (isempty(x))
    x = tt_rand(n,d,[K*ones(d,1);1],-1);
    rx = x.r;
    x = core2cell(x);
else
    if (x.d~=d)||(any(x.n~=n))
        error('Mismatching TT dimensionality of x0');
    end;
    rx = x.r;
    x = core2cell(x);
end;

% Parse auxiliary matrices, if any
if (~isempty(auxmat))
    raux = cell(K,1);
    maux = zeros(K,d); % Column sizes of aux
    for ki=1:K
        if (~isempty(auxmat{ki}))
            if (isa(auxmat{ki}, 'tt_matrix'))
                % Copy TT ranks and blocks from a tt_matrix
                raux{ki} = auxmat{ki}.r;
                maux(ki,:) = auxmat{ki}.m;
                if (auxmat{ki}.d~=d)||(any(auxmat{ki}.n~=n))
                    error('Mismatching TT dimensionalities in Aux matrix (%d)', ki);
                end;                
                auxmat{ki} = core2cell(auxmat{ki});
            else
                R = size(auxmat{ki},2);
                if (size(auxmat{ki},1)~=d)
                    error('Mismatching TT dimensionalities in Aux matrix (%d)', ki);
                end;
                raux{ki}=ones(d+1,R);
                for j=1:R
                    for i=1:d
                        if (issparse(auxmat{ki}{i,j}))
                            % Sparse matrix can only have TT ranks 1
                            [n1,n2]=size(auxmat{ki}{i,j});
                            ri = 1; ri2 = 1;
                        else
                            [ri,n1,n2,ri2]=size(auxmat{ki,kj}{i,j});
                        end;
                        if (n1~=n(i))
                            error('Mismatching TT dimensionalities in Aux matrix (%d)', ki);
                        end;
                        raux{ki}(i,j) = ri;
                        raux{ki}(i+1,j) = ri2;
                        maux(ki,i) = n2;
                    end;
                end;
            end;
        end;        
    end;
else
    auxmat = cell(K,1);
end;

% We need to normalize the components of X to compress them fairly.
% Store the norms here
scales = ones(K,1);

% Initialize the storage for residual Z
rz = [1; kickrank*ones(d-1,1); 1];
z = cell(d,1);

% An index showing which dimension introduces the RHS correction. No
% correction by default
dy_dim = 0;

% Reductions
XAX = cell(K,K); % Projections X_m'*A{i,j}*X_m of submatrices of A
XY = cell(K,1); % Similarly for y and z
ZAX = cell(K,K);
ZY = cell(K,1);
% Initialize the arrays
for kj=1:K
    for ki=1:K
        if (~isempty(A{ki,kj}))
            R = size(A{ki,kj}, 2);            
            XAX{ki,kj} = cell(d+1,R);
            ZAX{ki,kj} = cell(d+1,R);
            % Initialize border reductions by ones
            for j=1:R
                XAX{ki,kj}{1,j}=1; XAX{ki,kj}{d+1,j}=1;
                ZAX{ki,kj}{1,j}=1; ZAX{ki,kj}{d+1,j}=1;
            end;
        end;
    end;
    XY{kj} = cell(d+1,1);
    XY{kj}{1}=1; XY{kj}{d+1}=1;
    ZY{kj} = cell(d+1,1);
    ZY{kj}{1}=1; ZY{kj}{d+1}=1;
end;
% Reductions of aux matrices
XBX = cell(K,1);
for ki=1:K
    if (~isempty(auxmat{ki}))
        R = size(auxmat{ki},2);
        XBX{ki} = cell(d+1,R);
        for j=1:R
            XBX{ki}{1,j}=1; XBX{ki}{d+1,j}=1;
        end;
    end;
end;

% Orthogonalization. We have to check, where is the block enumerator
if (rx(d+1)>1)
    dir=-1;
else
    dir=1;
end;
if (dir>0)
    for i=d:-1:2
        % Orthogonalization itself, right to left
        [x{i-1}, x{i}, rx(i)]=orthogonalize_right(x{i-1}, x{i}, rx(i-1), n(i-1), rx(i), n(i), rx(i+1));
        % Create Z
        z{i} = randn(n(i)*rz(i+1), rz(i));
        [z{i},~]=qr(z{i},0);
        rz(i) = size(z{i},2);
        z{i} = z{i}.';
        % Reductions
        for kj=1:K
            for ki=1:K
                % Reductions with A
                if (~isempty(A{ki,kj}))
                    R = size(A{ki,kj},2); % Canonical rank
                    XAX{ki,kj}(i,:) = rightreduce_matrix(XAX{ki,kj}(i+1,:), x{i}, A{ki,kj}(i,:), x{i}, ...
                        rx(i),n(i),rx(i+1), R,ra{ki,kj}(i,:),ra{ki,kj}(i+1,:), rx(i),n(i),rx(i+1));
                    % the same reduction for Z
                    ZAX{ki,kj}(i,:) = rightreduce_matrix(ZAX{ki,kj}(i+1,:), z{i}, A{ki,kj}(i,:), x{i}, ...
                        rz(i),n(i),rz(i+1), R,ra{ki,kj}(i,:),ra{ki,kj}(i+1,:), rx(i),n(i),rx(i+1));
                end;
            end;
            % Reductions with Y.
            XY{kj}(i) = rightreduce_vector(XY{kj}(i+1), x{i}, y{kj}(i), ...
                rx(i),n(i),rx(i+1), 1,ry{kj}(i),ry{kj}(i+1));
            ZY{kj}(i) = rightreduce_vector(ZY{kj}(i+1), z{i}, y{kj}(i), ...
                rz(i),n(i),rz(i+1), 1,ry{kj}(i),ry{kj}(i+1));
            % Reductions with aux
            if (~isempty(auxmat{kj}))&&(n(i)==maux(kj,i))
                R = size(auxmat{kj},2); % Canonical rank
                XBX{kj}(i,:) = rightreduce_matrix(XBX{kj}(i+1,:), x{i}, auxmat{kj}(i,:), x{i}, ...
                    rx(i),n(i),rx(i+1), R,raux{kj}(i,:),raux{kj}(i+1,:), rx(i),n(i),rx(i+1));
            end;
        end;
    end;
    
    % Reshape the first TT core to make the blocking size the 4th index
    x{1} = permute(x{1}, [2,3,1]);
    x{1} = reshape(x{1}, 1, n(1), rx(2), K);
    rx(1) = 1;
    
    i = 1;
else
    for i=1:d-1
        % Orthogonalization itself, left to right
        [x{i}, x{i+1}, rx(i+1)]=orthogonalize_left(x{i}, x{i+1}, rx(i), n(i), rx(i+1), n(i+1), rx(i+2));
        % Create Z
        z{i} = randn(rz(i)*n(i), rz(i+1));
        [z{i},~]=qr(z{i},0);
        rz(i+1) = size(z{i},2);
        % Reductions
        for kj=1:K
            for ki=1:K
                % Reductions with A
                if (~isempty(A{ki,kj}))
                    R = size(A{ki,kj},2); % Canonical rank
                    XAX{ki,kj}(i+1,:) = leftreduce_matrix(XAX{ki,kj}(i,:), x{i}, A{ki,kj}(i,:), x{i}, ...
                        rx(i),n(i),rx(i+1), R,ra{ki,kj}(i,:),ra{ki,kj}(i+1,:), rx(i),n(i),rx(i+1));
                    % the same reduction for Z
                    ZAX{ki,kj}(i+1,:) = leftreduce_matrix(ZAX{ki,kj}(i,:), z{i}, A{ki,kj}(i,:), x{i}, ...
                        rz(i),n(i),rz(i+1), R,ra{ki,kj}(i,:),ra{ki,kj}(i+1,:), rx(i),n(i),rx(i+1));
                end;
            end;
            % Reductions with Y.
            XY{kj}(i+1) = leftreduce_vector(XY{kj}(i), x{i}, y{kj}(i), ...
                rx(i),n(i),rx(i+1), 1,ry{kj}(i),ry{kj}(i+1));
            ZY{kj}(i+1) = leftreduce_vector(ZY{kj}(i), z{i}, y{kj}(i), ...
                rz(i),n(i),rz(i+1), 1,ry{kj}(i),ry{kj}(i+1));
            % Reductions with aux
            if (~isempty(auxmat{kj}))&&(maux(kj,i)==n(i))
                R = size(auxmat{kj},2); % Canonical rank
                XBX{kj}(i+1,:) = leftreduce_matrix(XBX{kj}(i,:), x{i}, auxmat{kj}(i,:), x{i}, ...
                    rx(i),n(i),rx(i+1), R,raux{kj}(i,:),raux{kj}(i+1,:), rx(i),n(i),rx(i+1));
            end;
        end;
    end;
    
    % Reshape the last TT core to make the blocking size the 4th index
    x{d} = reshape(x{d}, rx(d), n(d), 1, K);
    rx(d+1) = 1;
    i = d;
end;


% Main cycle
swp = 1;
max_res = 0;
max_dx = 0;
max_z = 0;
while (swp<=nswp)
    % The sizes of the local problem
    locpos = rx(i)*n(i)*rx(i+1)*ones(K,1);
    locpos = cumsum([1;locpos]);
    
    % Extract matrix parts for the local problem
    XY1 = cell(K,1);
    XY2 = cell(K,1);
    Ai = cell(K,K);
    XAX1 = cell(K,K);
    XAX2 = cell(K,K);
    ZAX1 = cell(K,K);
    ZAX2 = cell(K,K);
    XBX1 = cell(K,1);
    auxi = cell(K,1);
    XBX2 = cell(K,1);
    raK = cell(K,K);
    yi = cell(K,1);
    for kj=1:K
        for ki=1:K
            if (~isempty(A{ki,kj}))
                for j=1:size(A{ki,kj},2)
                    Ai{ki,kj} = A{ki,kj}(i,:);
                    XAX1{ki,kj} = XAX{ki,kj}(i,:);
                    XAX2{ki,kj} = XAX{ki,kj}(i+1,:);
                    ZAX1{ki,kj} = ZAX{ki,kj}(i,:);
                    ZAX2{ki,kj} = ZAX{ki,kj}(i+1,:);
                    raK{ki,kj} = ra{ki,kj}([i,i+1],:);
                end;
            end;
        end;        
        if (i==dy_dim)
            % Working with the correcting block. Only the first term is the
            % "real" RHS
            XY1{kj} = XY{kj}(i,1);
            XY2{kj} = XY{kj}(i+1,1);
            yi{kj} = y{kj}(i,1);
        else
            % Otherwise, take all of them
            XY1{kj} = XY{kj}(i,:);
            XY2{kj} = XY{kj}(i+1,:);
            yi{kj} = y{kj}(i,:);
        end;
        if (~isempty(auxmat{kj}))
            XBX1{kj} = XBX{kj}(i,:);
            XBX2{kj} = XBX{kj}(i+1,:);
            auxi{kj} = auxmat{kj}(i,:);
        end;
    end;
    
    % Extract the RHS and the initial guess
    rhs_shf = zeros(locpos(K+1)-1, 1);
    for ki=1:K
        % Shifted RHS if there was correction -- for controlling the
        % residual
        R = size(y{ki},2);
        rhs_shf(locpos(ki):locpos(ki+1)-1) = assemble_local_vector(XY{ki}(i,:), y{ki}(i,:), XY{ki}(i+1,:), R,ry{ki}(i,:),ry{ki}(i+1,:), rx(i),n(i),rx(i+1));
    end;
    sol_prev = x{i}(:); % the size should be r1*n*r2*K already
    
    % Check the previous residual
    res_prev = norm(block_local_matvec(sol_prev, XAX1, Ai, XAX2, rx(i),n(i),rx(i+1), rx(i),n(i),rx(i+1), raK)-rhs_shf)/norm(rhs_shf);
    df = [];    
    if (res_prev>tol/sqrt(d)/resid_damp)||(i==dy_dim)
        % Solve the problem if needed
        if (all(cellfun('isempty', auxi)))
            sol = solfun(i, XAX1,Ai,XAX2, XY1,yi,XY2, tol/sqrt(d)/resid_damp, sol_prev);
        else
            sol = solfun(i, XAX1,Ai,XAX2, XY1,yi,XY2, tol/sqrt(d)/resid_damp, sol_prev, XBX1,auxi,XBX2);
        end;
        if (isa(sol,'cell'))
            % The RHS correction came
            df = sol{2};
            sol = sol{1};
        end;
        if (isempty(sol))
            % The user's procedure haven't solved the problem
            % Try the general approach
            sol = loc_solve_default(XAX1, Ai, XAX2, XY1, yi, XY2, tol/sqrt(d)/resid_damp, sol_prev, local_iters, max_full_size);
        end;
    else
        sol = sol_prev;
    end;
    
    % If the RHS correction came
    % format: df = [r x m x r] is the omitted solution
    if (~isempty(df))
        dy_dim = i;
        
        for ki=1:K
            % Populate y terms by -Bx
            R = size(auxmat{ki}, 2);
            for m=1:R
                for j=1:d
                    % Perform TT-matvec
                    if (j==i)
                        u = reshape(df, rx(j), maux(ki,j)*rx(j+1));
                        u = -u;
                    else
                        u = reshape(x{j}, rx(j), n(j)*rx(j+1));
                    end;
                    u = u.';
                    u = reshape(u, maux(ki,j), rx(j+1)*rx(j));
                    auxi = auxmat{ki}{j};
                    if (~issparse(auxi))
                        auxi = reshape(auxi, raux{ki}(j,m)*n(j)*maux(ki,j), raux{ki}(j+1,m));
                        auxi = auxi.';
                        auxi = reshape(auxi, raux{ki}(j+1,m)*raux{ki}(j,m)*n(j), maux(ki,j));
                    end;
                    u = auxi*u;
                    u = reshape(u, raux{ki}(j+1,m), raux{ki}(j,m), n(j), rx(j+1), rx(j));
                    u = permute(u, [2,5,3,1,4]);
                    u = reshape(u, raux{ki}(j,m)*rx(j), n(j), raux{ki}(j+1,m)*rx(j+1));
                    y{ki}{j,m+1} = u;
                    ry{ki}(j,m+1) = raux{ki}(j,m)*rx(j);
                    ry{ki}(j+1,m+1) = raux{ki}(j+1,m)*rx(j+1);
                end;
            end;
            % Rebuild reductions
            XY{ki}(1,1:R+1) = num2cell(ones(1,R+1));
            XY{ki}(d+1,1:R+1) = num2cell(ones(1,R+1));
            ZY{ki}(1,1:R+1) = num2cell(ones(1,R+1));
            ZY{ki}(d+1,1:R+1) = num2cell(ones(1,R+1));
            for j=1:i-1
                XY{ki}(j+1,2:R+1) = leftreduce_vector(XY{ki}(j,2:R+1), x{j}, y{ki}(j,2:R+1), ...
                    rx(j),n(j),rx(j+1), R,ry{ki}(j,2:R+1),ry{ki}(j+1,2:R+1));
                ZY{ki}(j+1,2:R+1) = leftreduce_vector(ZY{ki}(j,2:R+1), z{j}, y{ki}(j,2:R+1), ...
                    rz(j),n(j),rz(j+1), R,ry{ki}(j,2:R+1),ry{ki}(j+1,2:R+1));
            end;
            for j=d:-1:i+1
                XY{ki}(j,2:R+1) = rightreduce_vector(XY{ki}(j+1,2:R+1), x{j}, y{ki}(j,2:R+1), ...
                    rx(j),n(j),rx(j+1), R,ry{ki}(j,2:R+1),ry{ki}(j+1,2:R+1));
                ZY{ki}(j,2:R+1) = rightreduce_vector(ZY{ki}(j+1,2:R+1), z{j}, y{ki}(j,2:R+1), ...
                    rz(j),n(j),rz(j+1), R,ry{ki}(j,2:R+1),ry{ki}(j+1,2:R+1));
            end;
        end;
    end;
    
    % Normalize solutions
    sol = reshape(sol, rx(i)*n(i)*rx(i+1), K);
    for kj=1:K
        scales(kj) = norm(sol(:,kj));
        sol(:,kj) = sol(:,kj)/scales(kj);
    end;

    % Rebuild the shifted RHS, we need it to control the residual
    rhs_shf = zeros(locpos(K+1)-1, 1);
    for ki=1:K
        rhs_shf(locpos(ki):locpos(ki+1)-1) = assemble_local_vector(XY{ki}(i,:), y{ki}(i,:), XY{ki}(i+1,:), size(y{ki},2),ry{ki}(i,:),ry{ki}(i+1,:), rx(i),n(i),rx(i+1));
    end;
    res_new = norm(block_local_matvec(sol*diag(scales), XAX1, Ai, XAX2, rx(i),n(i),rx(i+1), rx(i),n(i),rx(i+1), raK)-rhs_shf)/norm(rhs_shf);
    
    fprintf('--- swp=%d, i=%d, dx=[ ', swp, i);
    % Measure the error
    sol_prev = reshape(sol_prev, rx(i)*n(i)*rx(i+1), K);
    for kj=1:K
        dx = norm(sol(:,kj)-sol_prev(:,kj)/scales(kj));
        max_dx = max(max_dx, dx);
        fprintf('%3.3e ', dx);
    end;
    fprintf('], res_prev=%3.3e, res_new=%3.3e\n', res_prev, res_new);
    max_res = max(max_res, res_prev);
    
    % Save old x ranks. We will need them to compute the residual
    rx1old = rx(i);
    rx2old = rx(i+1);
    % Truncation and enrichment -- dependent on the direction
    if (dir>0)&&(i<d)
        % Truncation
        sol = reshape(sol, rx(i)*n(i), rx(i+1)*K);
        [u,s,v]=svd(sol, 'econ');
        s = diag(s);
        % Residual thresholding
        sol_prev = u(:,1)*s(1)*v(:,1)';
        sol_prev = reshape(sol_prev, rx(i)*n(i)*rx(i+1), K);
        sol_prev = sol_prev*diag(scales);
        resid = norm(block_local_matvec(sol_prev, XAX1, Ai, XAX2, rx(i),n(i),rx(i+1), rx(i),n(i),rx(i+1), raK)-rhs_shf)/norm(rhs_shf);
        for rnew=2:numel(s)
            if (resid>max(tol/sqrt(d), res_new*resid_damp))
                sol_new_solid = u(:,rnew)*s(rnew)*v(:,rnew)';
                sol_new_solid = reshape(sol_new_solid, rx(i)*n(i)*rx(i+1), K);
                sol_prev = sol_prev+sol_new_solid*diag(scales);
                resid = norm(block_local_matvec(sol_prev, XAX1, Ai, XAX2, rx(i),n(i),rx(i+1), rx(i),n(i),rx(i+1), raK)-rhs_shf)/norm(rhs_shf);
            else
                break;
            end;
        end;
        u = u(:,1:rnew);
        s = s(1:rnew);
        v = v(:,1:rnew);
        v = diag(s)*v';
        
        % Prepare the solid truncated solution -- to compute  Z
        sol_new_solid = u*v;
        sol_new_solid = reshape(sol_new_solid, rx(i)*n(i)*rx(i+1), K);
        sol_new_solid = sol_new_solid*diag(scales);
        
        % Compute the enrichment
        % Solid Z
        z_Ax = block_local_matvec(sol_new_solid, XAX1, Ai, ZAX2, rx(i),n(i),rz(i+1), rx(i),n(i),rx(i+1), raK);
        locposZ = rx(i)*n(i)*rz(i+1)*ones(K,1);
        locposZ = cumsum([1;locposZ]);
        z_y = zeros(locposZ(K+1)-1, 1);
        for ki=1:K
            z_y(locposZ(ki):locposZ(ki+1)-1) = assemble_local_vector(XY{ki}(i,:), y{ki}(i,:), ZY{ki}(i+1,:), size(y{ki},2),ry{ki}(i,:),ry{ki}(i+1,:), rx(i),n(i),rz(i+1));
        end;
        % Concatenate the residuals
        z_loc = zeros(rx(i)*n(i)*rz(i+1), K);
        for ki=1:K
            z_new = z_y(locposZ(ki):locposZ(ki+1)-1) - z_Ax(locposZ(ki):locposZ(ki+1)-1);
            if (norm(z_new, 'fro')>0)
                z_new = z_new/norm(z_new, 'fro');
            end;
            z_loc(:,ki) = z_new;
        end;
        % Take the SVD of the residual
        z_loc = reshape(z_loc, rx(i)*n(i), rz(i+1)*K);
        [z_loc,~,~]=svd(z_loc,'econ');
        z_loc = z_loc(:,1:min(kickrank,size(z_loc,2)));
        % Perform the enrichment and QR it
        [u,rv]=qr([u, z_loc], 0);
        rv = rv(:,1:size(v,1));
        v = rv*v; % its size is r2new x r2
        rnew = size(u, 2);
        % Prepare the next TT block. Steven White loves us =)
        v = reshape(v, rnew*rx(i+1), K);
        % don't forget the norms!
        v = v*diag(scales);
        v = v.';
        v = reshape(v, K*rnew, rx(i+1));
        v = v*reshape(x{i+1}, rx(i+1), n(i+1)*rx(i+2));
        v = reshape(v, K, rnew*n(i+1)*rx(i+2));
        v = v.';
        % Update the rank
        rx(i+1) = rnew;
        % Put back to X
        x{i} = reshape(u, rx(i), n(i), rx(i+1));
        x{i+1} = reshape(v, rx(i+1), n(i+1), rx(i+2), K);
        
        % Compute the reductions of X
        for kj=1:K
            for ki=1:K
                % ... with A
                if (~isempty(A{ki,kj}))
                    R = size(A{ki,kj},2); % Canonical rank
                    XAX{ki,kj}(i+1,:) = leftreduce_matrix(XAX{ki,kj}(i,:), x{i}, Ai{ki,kj}, x{i}, ...
                        rx(i),n(i),rx(i+1), R,ra{ki,kj}(i,:),ra{ki,kj}(i+1,:), rx(i),n(i),rx(i+1));
                end;
            end;
            % Reductions with Y
            XY{kj}(i+1,:) = leftreduce_vector(XY{kj}(i,:), x{i}, y{kj}(i,:), ...
                rx(i),n(i),rx(i+1), size(y{kj},2),ry{kj}(i,:),ry{kj}(i+1,:));
            % Reductions with aux
            if (~isempty(auxmat{kj}))
                R = size(auxmat{kj},2); % Canonical rank
                if (maux(kj,i)==n(i))&&(~isempty(XBX{kj}{i,1}))
                    XBX{kj}(i+1,:) = leftreduce_matrix(XBX{kj}(i,:), x{i}, auxmat{kj}(i,:), x{i}, ...
                        rx(i),n(i),rx(i+1), R,raux{kj}(i,:),raux{kj}(i+1,:), rx(i),maux(kj,i),rx(i+1));
                else
                    % matrix is not square -> reduction is not possible
                    XBX{kj}(i+1,:) = cell(1,R);
                end;
            end;
        end;  
        
    elseif (dir<0)&&(i>1)
        % Truncation
        sol = sol.';
        sol = reshape(sol, K*rx(i), n(i)*rx(i+1));
        [u,s,v]=svd(sol, 'econ');
        s = diag(s);
        % Residual thresholding
        sol_prev = u(:,1)*s(1)*v(:,1)';
        sol_prev = reshape(sol_prev, K, rx(i)*n(i)*rx(i+1));
        sol_prev = sol_prev.'*diag(scales);
        resid = norm(block_local_matvec(sol_prev, XAX1, Ai, XAX2, rx(i),n(i),rx(i+1), rx(i),n(i),rx(i+1), raK)-rhs_shf)/norm(rhs_shf);
        for rnew=2:numel(s)
            if (resid>max(tol/sqrt(d), res_new*resid_damp))
                sol_new_solid = u(:,rnew)*s(rnew)*v(:,rnew)';
                sol_new_solid = reshape(sol_new_solid, K, rx(i)*n(i)*rx(i+1));
                sol_new_solid = sol_new_solid.';
                sol_prev = sol_prev+sol_new_solid*diag(scales);
                resid = norm(block_local_matvec(sol_prev, XAX1, Ai, XAX2, rx(i),n(i),rx(i+1), rx(i),n(i),rx(i+1), raK)-rhs_shf)/norm(rhs_shf);
            else
                break;
            end;
        end;
        u = u(:,1:rnew);
        s = s(1:rnew);
        v = v(:,1:rnew);
        v = v';
        u = u*diag(s);
        
        % Prepare the solid truncated solution -- to compute Z
        sol_new_solid = u*v;
        sol_new_solid = reshape(sol_new_solid, K, rx(i)*n(i)*rx(i+1));
        sol_new_solid = sol_new_solid.'*diag(scales);
        
        % Compute the enrichment
        % Solid Z
        z_Ax = block_local_matvec(sol_new_solid, ZAX1, Ai, XAX2, rz(i),n(i),rx(i+1), rx(i),n(i),rx(i+1), raK);
        locposZ = rz(i).*n(i).*rx(i+1)*ones(K,1);
        locposZ = cumsum([1;locposZ]);
        z_y = zeros(locposZ(K+1)-1, 1);
        for ki=1:K
            z_y(locposZ(ki):locposZ(ki+1)-1) = assemble_local_vector(ZY{ki}(i,:), y{ki}(i,:), XY{ki}(i+1,:), size(y{ki},2),ry{ki}(i,:),ry{ki}(i+1,:), rz(i),n(i),rx(i+1));
        end;
        % Concatenate the residuals
        z_loc = zeros(K, rz(i)*n(i)*rx(i+1));
        for ki=1:K
            z_new = z_y(locposZ(ki):locposZ(ki+1)-1) - z_Ax(locposZ(ki):locposZ(ki+1)-1);
            if (norm(z_new, 'fro')>0)
                z_new = z_new/norm(z_new, 'fro');
            end;
            z_loc(ki,:) = z_new.';
        end;
        % Take the SVD of the residual
        z_loc = reshape(z_loc, K*rz(i), n(i)*rx(i+1));
        [~,~,z_loc]=svd(z_loc,'econ');
        z_loc = z_loc(:,1:min(kickrank,size(z_loc,2)))';
        % Perform the enrichment and QR it
        [v,rv]=qr([v; z_loc].', 0);
        rv = rv(:,1:size(u,2));
        u = u*rv.'; % its size is r1 x r1new
        v = v.';
        rnew = size(v, 1);
        % Prepare the next TT block. Steven White loves us =)
        u = reshape(u, K, rx(i)*rnew);
        % don't forget the norms!
        u = diag(scales)*u;
        u = u.';
        u = reshape(u, rx(i), rnew*K);
        u = reshape(x{i-1}, rx(i-1)*n(i-1), rx(i))*u;
        u = reshape(u, rx(i-1), n(i-1), rnew, K);
        % Update the rank
        rx(i) = rnew;
        % Put back to X
        x{i} = reshape(v, rx(i), n(i), rx(i+1));
        x{i-1} = u;
        
        % Compute the reductions of X
        for kj=1:K
            % with A
            for ki=1:K
                if (~isempty(A{ki,kj}))
                    R = size(A{ki,kj},2); % Canonical rank
                    XAX{ki,kj}(i,:) = rightreduce_matrix(XAX{ki,kj}(i+1,:), x{i}, Ai{ki,kj}, x{i}, ...
                        rx(i),n(i),rx(i+1), R,ra{ki,kj}(i,:),ra{ki,kj}(i+1,:), rx(i),n(i),rx(i+1));
                end;
            end;
            % Reductions with Y
            XY{kj}(i,:) = rightreduce_vector(XY{kj}(i+1,:), x{i}, y{kj}(i,:), ...
                rx(i),n(i),rx(i+1), size(y{kj},2),ry{kj}(i,:),ry{kj}(i+1,:));
            % Reductions with aux
            if (~isempty(auxmat{kj}))
                R = size(auxmat{kj},2); % Canonical rank
                if (maux(kj,i)==n(i))&&(~isempty(XBX{kj}{i+1,1}))
                    XBX{kj}(i,:) = rightreduce_matrix(XBX{kj}(i+1,:), x{i}, auxmat{kj}(i,:), x{i}, ...
                        rx(i),n(i),rx(i+1), R,raux{kj}(i,:),raux{kj}(i+1,:), rx(i),maux(kj,i),rx(i+1));
                else
                    XBX{kj}(i,:) = cell(1,R);
                end;
            end;            
        end;
        
    else % border TT block
        % Just copy sol to the TT format storage
        sol = sol*diag(scales);
        x{i} = reshape(sol, rx(i), n(i), rx(i+1), K);
        sol_new_solid = sol;
    end;
    
    % Update the residual
    % Solid Z
    z_Ax = block_local_matvec(sol_new_solid, ZAX1, Ai, ZAX2, rz(i),n(i),rz(i+1), rx1old,n(i),rx2old, raK);
    locposZ = rz(i).*n(i).*rz(i+1)*ones(K,1);
    locposZ = cumsum([1;locposZ]);
    z_y = zeros(locposZ(K+1)-1, 1);
    for ki=1:K
        z_y(locposZ(ki):locposZ(ki+1)-1) = assemble_local_vector(ZY{ki}(i,:), y{ki}(i,:), ZY{ki}(i+1,:), size(y{ki},2),ry{ki}(i,:),ry{ki}(i+1,:), rz(i),n(i),rz(i+1));
    end;
    % Cycle over residual blocks
    z{i} = [];
    for ki=1:K
        z_new = z_y(locposZ(ki):locposZ(ki+1)-1) - z_Ax(locposZ(ki):locposZ(ki+1)-1);
        max_z = max(max_z, norm(z_new));
        if (norm(z_new, 'fro')>0)
            z_new = z_new/norm(z_new, 'fro');
        end;
        % Truncate the residual
        if (dir>0)
            z_new = reshape(z_new, rz(i)*n(i), rz(i+1));
        else
            z_new = reshape(z_new, rz(i), n(i)*rz(i+1));
            z_new = z_new.';
        end;
        z{i} = [z{i}, z_new];
    end;
    
    % Orthogonalize the residuals
    [z{i},~,~]=svd(z{i},'econ');
    z{i} = z{i}(:, 1:min(kickrank,size(z{i},2)));
    if (dir>0)&&(i<d)
        rz(i+1) = size(z{i},2);
    elseif (dir<0)&&(i>1)
        rz(i) = size(z{i},2);
        z{i} = z{i}.';
    end;
    
    % Compute reductions of Z
    if (dir>0)&&(i<d)
        for kj=1:K
            % with A
            for ki=1:K
                if (~isempty(A{ki,kj}))
                    R = size(A{ki,kj},2); % Canonical rank
                    ZAX{ki,kj}(i+1,:) = leftreduce_matrix(ZAX{ki,kj}(i,:), z{i}, Ai{ki,kj}, x{i}, ...
                        rz(i),n(i),rz(i+1), R,ra{ki,kj}(i,:),ra{ki,kj}(i+1,:), rx(i),n(i),rx(i+1));
                end;
            end;
            % Reductions with Y
            ZY{kj}(i+1,:) = leftreduce_vector(ZY{kj}(i,:), z{i}, y{kj}(i,:), ...
                rz(i),n(i),rz(i+1), size(y{kj},2),ry{kj}(i,:),ry{kj}(i+1,:));
        end;
    elseif (dir<0)&&(i>1)
        for kj=1:K
            % with A
            for ki=1:K
                if (~isempty(A{ki,kj}))
                    R = size(A{ki,kj},2); % Canonical rank
                    ZAX{ki,kj}(i,:) = rightreduce_matrix(ZAX{ki,kj}(i+1,:), z{i}, Ai{ki,kj}, x{i}, ...
                        rz(i),n(i),rz(i+1), R,ra{ki,kj}(i,:),ra{ki,kj}(i+1,:), rx(i),n(i),rx(i+1));
                end;
            end;
            % Reductions with Y
            ZY{kj}(i,:) = rightreduce_vector(ZY{kj}(i+1,:), z{i}, y{kj}(i,:), ...
                rz(i),n(i),rz(i+1), size(y{kj},2),ry{kj}(i,:),ry{kj}(i+1,:));
        end;
    end;
    
    
    % Reversing, end of sweep
    i = i+dir;
    if (i>d)||(i<1)
        
        if ((exitdir>0)&&(i>d))||((exitdir<0)&&(i<1))
            % Full backward-forward sweep done
            fprintf('==amen_block_solve== swp=%d, max_dx=%3.3e, max_res=%3.3e, max_rank=%d, max_z=%3.3e\n', swp, max_dx, max_res, max(rx(:)), max_z);
            if (max_res<tol_exit)
                break;
            end;
            
            max_dx = 0;
            max_res = 0;
            max_z = 0;
            swp = swp+1;
        else
            % Only half-sweep is finished
            fprintf('\t swp=%d.5, max_dx=%3.3e, max_res=%3.3e, max_rank=%d, max_z=%3.3e\n', swp-1, max_dx, max_res, max(rx(:)), max_z);
        end;
        
        dir = -dir;
        i = i+dir;
    end;
end;

% Reshape the last block to the BTT format
if (exitdir>0)
    x{d} = reshape(x{d}, rx(d), n(d), K);
else
    x{1} = reshape(x{1}, n(1)*rx(2), K);
    x{1} = x{1}.';
    x{1} = reshape(x{1}, K, n(1), rx(2));
end;
x = cell2core(tt_tensor, x);
end


% R-L Orthogonalization, from x1 to x0
function [x0,x1,r1]=orthogonalize_right(x0,x1,r0,n0,r1,n1,r2)
x1 = reshape(x1, r1, n1*r2);
[x1,rv]=qr(x1.', 0);
x0 = reshape(x0, r0*n0, r1);
x0 = x0*rv.';
r1 = size(x1,2);
x1 = x1.';
x1 = reshape(x1, r1, n1, r2);
x0 = reshape(x0, r0, n0, r1);
end

% L-R Orthogonalization, from x1 to x2
function [x1,x2,r2]=orthogonalize_left(x1,x2,r1,n1,r2,n2,r3)
x1 = reshape(x1, r1*n1, r2);
[x1,rv]=qr(x1, 0);
x2 = reshape(x2, r2, n2*r3);
x2 = rv*x2;
r2 = size(x1,2);
x1 = reshape(x1, r1, n1, r2);
x2 = reshape(x2, r2, n2, r3);
end

% Accumulates the left reduction W{1:k}'*A{1:k}*X{1:k}
function [WAX2] = leftreduce_matrix(WAX1, w, A, x, rw1,n,rw2, Ra,ra1,ra2, rx1,m,rx2)
% Left WAX has the form of the first matrix TT block, i.e. [rw, rx, ra]
WAX2 = WAX1;
wc = reshape(w, rw1, n*rw2);
xc = reshape(x, rx1*m, rx2);
for k=1:Ra
    WAX2{k} = reshape(WAX2{k}, rw1, rx1*ra1(k));
    WAX2{k} = wc'*WAX2{k}; % size n rw2 x rx1 ra1
    WAX2{k} = reshape(WAX2{k}, n, rw2*rx1*ra1(k));
    WAX2{k} = WAX2{k}.';
    WAX2{k} = reshape(WAX2{k}, rw2*rx1, ra1(k)*n);
    tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
    WAX2{k} = WAX2{k}*tmp; % size rw2 rx1 m ra2
    WAX2{k} = reshape(WAX2{k}, rw2, rx1*m*ra2(k));
    WAX2{k} = WAX2{k}.';
    WAX2{k} = reshape(WAX2{k}, rx1*m, ra2(k)*rw2);
    WAX2{k} = xc.'*WAX2{k}; % size rx2, ra2 rw2
    WAX2{k} = reshape(WAX2{k}, rx2*ra2(k), rw2);
    WAX2{k} = WAX2{k}.';
    WAX2{k} = reshape(WAX2{k}, rw2, rx2, ra2(k));
end;
end

% Accumulates the left reduction W{1:k}'*X{1:k}
function [WX2] = leftreduce_vector(WX1, w, x, rw1,n,rw2, Rx,rx1,rx2)
% Left WX has the form of the first vector TT block, i.e. [rw, rx]
WX2 = WX1;
wc = reshape(w, rw1, n*rw2);
for k=1:Rx
    WX2{k} = wc'*WX2{k}; % size n rw2 x rx1
    WX2{k} = reshape(WX2{k}, n, rw2*rx1(k));
    WX2{k} = WX2{k}.';
    WX2{k} = reshape(WX2{k}, rw2, rx1(k)*n);
    tmp = reshape(x{k}, rx1(k)*n, rx2(k));
    WX2{k} = WX2{k}*tmp; % size rw2, rx2
end;
end

% Accumulates the right reduction W{k:d}'*A{k:d}*X{k:d}
function [WAX1] = rightreduce_matrix(WAX2, w, A, x, rw1,n,rw2, Ra,ra1,ra2, rx1,m,rx2)
% Right WAX has the form of the last matrix TT block, i.e. [ra, rw, rx]
WAX1 = WAX2;
wc = reshape(w, rw1, n*rw2);
wc = conj(wc);
xc = reshape(x, rx1*m, rx2);
for k=1:Ra
    WAX1{k} = reshape(WAX1{k}, ra2(k)*rw2, rx2);
    WAX1{k} = xc*WAX1{k}.'; % size rx1 m x ra2 rw2
    WAX1{k} = reshape(WAX1{k}, rx1, m*ra2(k)*rw2);
    WAX1{k} = WAX1{k}.';
    WAX1{k} = reshape(WAX1{k}, m*ra2(k), rw2*rx1);
    tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
    WAX1{k} = tmp*WAX1{k}; % size ra1(k)*n, rw2*rx1
    WAX1{k} = reshape(WAX1{k}, ra1(k), n*rw2*rx1);
    WAX1{k} = WAX1{k}.';
    WAX1{k} = reshape(WAX1{k}, n*rw2, rx1*ra1(k));
    WAX1{k} = wc*WAX1{k}; % size rw1, rx1 ra1
    WAX1{k} = reshape(WAX1{k}, rw1*rx1, ra1(k));
    WAX1{k} = WAX1{k}.';
    WAX1{k} = reshape(WAX1{k}, ra1(k), rw1, rx1);
end;
end


% Accumulates the right reduction W{k:d}'*X{k:d}
function [WX1] = rightreduce_vector(WX2, w, x, rw1,n,rw2, Rx,rx1,rx2)
% Right WX has the form of the last vector TT block, i.e. [rx, rw]
WX1 = WX2;
wc = reshape(w, rw1, n*rw2);
for k=1:Rx
    tmp = reshape(x{k}, rx1(k)*n, rx2(k));
    WX1{k} = tmp*WX1{k}; % size rx1 n x rw2
    WX1{k} = reshape(WX1{k}, rx1(k), n*rw2);
    WX1{k} = WX1{k}*wc'; % size rx1, rw1
end;
end

% A matrix-vectors product for the matrix in the 3D TT (WAX1-A-WAX2), and
% full vectors of size (rx1*m*rx2) x b. Returns (rw1*n*rw2) x b
function [w]=local_matvec(x, rx1,m,rx2,b, rw1,n,rw2, WAX1, A, WAX2, Ra,ra1,ra2)
w = zeros(rw1*n*rw2, b);
xc = reshape(x, rx1*m*rx2, b);
xc = xc.';
xc = reshape(xc, b*rx1*m, rx2);
for k=1:Ra
    tmp = reshape(WAX2{k}, ra2(k)*rw2, rx2);
    wk = xc*tmp.';
    wk = reshape(wk, b*rx1, m*ra2(k)*rw2);
    wk = wk.';
    wk = reshape(wk, m*ra2(k), rw2*b*rx1);
    tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
    wk = tmp*wk;
    wk = reshape(wk, ra1(k)*n*rw2*b, rx1);
    wk = wk.';
    wk = reshape(wk, rx1*ra1(k), n*rw2*b);
    tmp = reshape(WAX1{k}, rw1, rx1*ra1(k));
    wk = tmp*wk;
    wk = reshape(wk, rw1*n*rw2, b);
    w = w+wk;
end;
end

% Builds the full (rw1*n*rw2) x 1 vector from its TT blocks
function [w]=assemble_local_vector(WX1, x, WX2, Rx,rx1,rx2, rw1,n,rw2)
w = zeros(rw1*n*rw2, 1);
for k=1:Rx
    wk = reshape(x{k}, rx1(k), n*rx2(k));
    wk = WX1{k}*wk;
    wk = reshape(wk, rw1*n, rx2(k));
    wk = wk*WX2{k};
    wk = reshape(wk, rw1*n*rw2, 1);
    w = w+wk;
end;
end

% Block local Matvec, where y{i} = sum_j local_matvec(A{i,j}, x{j})
function [y] = block_local_matvec(x, WAX1, Ai, WAX2, rw1,n,rw2, rx1,m,rx2, ra)
K = size(Ai,1);
locposX = rx1*m*rx2*ones(K,1);
locposX = cumsum([1; locposX]);
locposY = rw1*n*rw2*ones(K,1);
locposY = cumsum([1; locposY]);
y = zeros(rw1*n*rw2*K,1);
for j=1:K
    xj = x(locposX(j):locposX(j+1)-1);
    for i=1:K
        if (~isempty(Ai{i,j}))
            R = size(Ai{i,j},2);
            yi = local_matvec(xj, rx1,m,rx2,1, rw1,n,rw2, WAX1{i,j}, Ai{i,j}, WAX2{i,j}, R,ra{i,j}(1,:),ra{i,j}(2,:));
            y(locposY(i):locposY(i+1)-1) = y(locposY(i):locposY(i+1)-1) + yi;
        end;
    end;
end;
end

function [y]=loc_solve_default(XAX1, Ai, XAX2, XY1,yi,XY2, tol, sol_prev, local_iters, max_full_size)
% Extract the sizes
K = size(Ai,1);
ra = cell(K,K);
for ki=1:K
    for kj=1:K
        if (~isempty(Ai{ki,kj}))
            R = size(Ai{ki,kj},2);
            r1 = size(XAX1{ki,kj}{1},1);
            r2 = size(XAX2{ki,kj}{1},2);
            ra{ki,kj} = ones(2,R);
            for j=1:R
                if (issparse(Ai{ki,kj}{j}))
                    ra{ki,kj}(:,j) = 1;
                    n = size(Ai{ki,kj}{1},1);
                else
                    ra{ki,kj}(1,j) = size(Ai{ki,kj}{j},1);
                    ra{ki,kj}(2,j) = size(Ai{ki,kj}{j},4);
                    n = size(Ai{ki,kj}{1},2);
                end;
            end;
        end;
    end;
end;

% Assemble the RHS
sz = r1*n*r2*ones(K,1);
locpos = cumsum([1; sz]);
rhs = zeros(r1*n*r2*K, 1);
for ki=1:K
    ry1 = cellfun(@(x)size(x,1), yi{ki});
    ry2 = cellfun(@(x)size(x,3), yi{ki});
    rhs(locpos(ki):locpos(ki+1)-1) = assemble_local_vector(XY1{ki}, yi{ki}, XY2{ki}, size(yi{ki},2),ry1,ry2, r1,n,r2);
end;

if (numel(rhs)<max_full_size)
    [B,sparseflag]=assemble_block_matrix(XAX1, Ai, XAX2, r1,n,r2, r1,n,r2, ra);
    y = rhs;
    if (sparseflag)
        % We'll put the mode to the senior index
        K = size(sz, 1);
        for ki=1:K
            yi = reshape(y(locpos(ki):locpos(ki+1)-1), r1*n, r2);
            yi = yi.';
            yi = reshape(yi, r2*r1*n, 1);
            y(locpos(ki):locpos(ki+1)-1) = yi;
        end;
    end;
    y = B\y;
    if (sparseflag)
        % Recover the order back
        for ki=1:K
            yi = reshape(y(locpos(ki):locpos(ki+1)-1), r2, r1*n);
            yi = yi.';
            yi = reshape(yi, r1*n*r2, 1);
            y(locpos(ki):locpos(ki+1)-1) = yi;
        end;
    end;
else
    [y,~] = gmres(@(x)block_local_matvec(x, XAX1, Ai, XAX2, r1,n,r2, r1,n,r2, ra), rhs, [], tol, local_iters, [], [], sol_prev);
end;
end





function [B,sparseflag]=assemble_local_matrix(WAX1, A, WAX2, Ra,ra1,ra2, rw1,n,rw2, rx1,m,rx2)
% Check the sparsity of the matrix blocks
sparseflag = true;
for k=1:Ra
    if (~issparse(A{k}))
        sparseflag=false;
    end;
end;
if (sparseflag)
    B = sparse(rw2*rw1*n, rx2*rx1*m); % reverse order !!!
    % Currently only canonical sparse matrices are allowed
    for k=1:Ra
        tmp = reshape(WAX2{k}, rw2, rx2);
        tmp = sparse(tmp);
        Bk = reshape(WAX1{k}, rw1, rx1);
        Bk = sparse(Bk);
        Bk = kron(Bk, tmp); % mind endiannes
        Bk = kron(A{k}, Bk); % mind endiannes
        B = B+Bk;
    end;
else
    B = zeros(rw1*n*rw2, rx1*m*rx2);
    for k=1:Ra
        Bk = reshape(WAX1{k}, rw1*rx1, ra1(k));
        tmp = reshape(A{k}, ra1(k), n*m*ra2(k));
        if (issparse(tmp))
            tmp = full(tmp);
        end;
        Bk = Bk*tmp;
        Bk = reshape(Bk, rw1, rx1, n, m*ra2(k));
        Bk = permute(Bk, [1,3,2,4]);
        Bk = reshape(Bk, rw1*n*rx1*m, ra2(k));
        tmp = reshape(WAX2{k}, ra2(k), rw2*rx2);
        Bk = Bk*tmp;
        Bk = reshape(Bk, rw1*n, rx1*m, rw2, rx2);
        Bk = permute(Bk, [1,3,2,4]);
        Bk = reshape(Bk, rw1*n*rw2, rx1*m*rx2);
        B = B+Bk;
    end;
end;
end

function [B,sparseflag]=assemble_block_matrix(WAX1, A, WAX2, rw1,n,rw2, rx1,m,rx2, raK)
K = size(A,1);
szw = rw1*n*rw2;
szx = rx1*m*rx2;

sparseflag = true;
B = [];
for kj=1:K
    Bj = [];
    for ki=1:K
        if (~isempty(A{ki,kj}))
            [Bij,sfij] = assemble_local_matrix(WAX1{ki,kj}, A{ki,kj}, WAX2{ki,kj}, size(A{ki,kj},2),raK{ki,kj}(1,:),raK{ki,kj}(2,:), rw1,n,rw2, rx1,m,rx2);
            Bj = [Bj; Bij];
            sparseflag = sparseflag & sfij;
        else
            if (sparseflag)
                Bj = [Bj; sparse(szw, szx)];
            else
                Bj = [Bj; zeros(szw, szx)];
            end;
        end;
    end;
    B = [B, Bj];
end;

end
