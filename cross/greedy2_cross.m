function [y,Jyl,Jyr,ilocl,ilocr]=greedy2_cross(n, fun, tol, varargin)
% Two-site greedy cross interpolation scheme.
%   [y,Jyl,Jyr,ilocl,ilocr]=greedy2_cross(n, fun, tol, varargin)
% Tries to interpolate the tensor with mode sizes n specified by the
% function fun up to the accuracy tol using the "basic" greedy restricted 
% cross interpolation algorithm in the TT format.
% The method sequentially adds one pivot to each superblock in a forward
% DMRG-type half-sweep (k=1,...,d).
% 
% The input n should be a vector of mode sizes of length d,
% fun = @(ind)fun(ind) is a sought function of index ind.
% By default, ind is an array of d indices, and the function fun should
% return one value of the corresponding tensor entry.
% To speed up the computations, set the optional parameter 'vec' to true,
% and provide the function which takes ind as an array of sizes M x d, and
% returns an array of M values.
%
% In addition to the indexwise function, one may provide the value-wise
% function of another tt_tensor (funcrs style) via optional parameters 
% 'aux', 'auxfun' (see below).
% The resulting tensor is the sum    y(ind)=fun(ind) + auxfun(aux(ind)).
%
% The first output parameter is the computed tt_tensor, the rest four
% parameters return the pivot indices (left-global, right-global, left-local, right-local) for
% testing purposes.
%
% Optional arguments are provided in the form
% 'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so on. 
% The list of option names and default values:
%       o nswp - maximal number of sweeps [20]
%       o tol_exit - stopping difference between consecutive iterations [tol]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o vec - whether fun can accept and return vectorized values [false]
%       o aux - a tt_tensor for auxiliary funcrs contribution []
%       o auxfun - an auxiliary function defined pointwise at the elements
%           of aux []
%       o 'locsearch' - an algorithm for the error pivoting in superblocks:
%           try a lottery of n*r random elements ('lot'), or
%           conduct two-dimensional maxvol ALS iteration ('als').
%           The first one uses less evaluations, but may be costly due to a
%           nonvectorized MATLAB loop ['lot']
%
% 
% This procedure implements Algorithm 2 from 
%   D. Savostyanov, http://arxiv.org/pdf/1305.1818v2.pdf.
% Development of the MATLAB version: S. Dolgov.
% Please send feedback to: {dmitry.savostyanov,sergey.v.dolgov}@gmail.com
%
%---------------------------

if (~isempty(varargin))
    v1 = varargin{1};
    if (isa(v1, 'cell'))
        varargin=v1;
    end;
end;
vars = varargin;

nswp = 20;
tol_exit = tol;
verb = 1;
vec = false;
aux = []; % Extra tt_tensor to pass into this cross
auxfun = []; % the total function equals fun(ind)+auxfun(aux(ind)), since simple aux(ind) via tt_tensor/subsref suxx.
% locsearch = 'als';
locsearch = 'lot';

i = 1;
while (i<length(vars))
    switch lower(vars{i})
        case 'nswp'
            nswp=vars{i+1};
        case 'tol_exit'
            tol_exit=vars{i+1};  
        case 'verb'
            verb = vars{i+1};      
        case 'vec'
            vec = vars{i+1};   
        case 'aux'
            aux = vars{i+1}; 
        case 'auxfun'
            auxfun = vars{i+1};             
        case 'locsearch'
            locsearch = vars{i+1};               
        otherwise
            warning('Option %s was not recognized', vars{i});
    end;
    i=i+2;
end;

d = numel(n);
ry = ones(d+1,1);
y = cell(d,1);

if (~isempty(aux))
    raux = aux.r;
    aux = core2cell(aux);
    phiauxl = cell(d+1,1); phiauxl{1}=1; phiauxl{d+1}=1;
    phiauxr = cell(d+1,1); phiauxr{1}=1; phiauxr{d+1}=1;
end;

% Factorized inverse interpolation matrix -- in the form U^{-1}, L^{-1}
mid_inv = cell(d+1,2); mid_inv{1,1}=1; mid_inv{1,2}=1; mid_inv{d+1,1}=1; mid_inv{d+1,2}=1;


Jyr = cell(d+1,1);
Jyl = cell(d+1,1);
ilocr = cell(d+1,1);
ilocl = cell(d+1,1);

evalcnt = 0;

% Start with some rand indices
for i=1:d-1
    ilocl{i+1} = rand(1,1)*ry(i)*n(i);
    ilocl{i+1} = round(ilocl{i+1});
    if (ilocl{i+1}==0)
        ilocl{i+1} = 1;
    end;
    Jyl{i+1} = indexmerge(Jyl{i}, (1:n(i))');
    Jyl{i+1} = Jyl{i+1}(ilocl{i+1},:);
end;
% Perform one sweep back-forth to fix the initial index
for i=d:-1:1
    J = indexmerge(Jyl{i}, (1:n(i))', Jyr{i+1});
    evalcnt = evalcnt + size(J,1);
    cry1 = autovecfun(fun, J, vec);
    if (~isempty(aux))
        craux1 = reshape(aux{i}, raux(i), n(i)*raux(i+1));
        craux1 = phiauxl{i}*craux1;
        craux1 = reshape(craux1, ry(i)*n(i), raux(i+1));
        craux1 = craux1*phiauxr{i+1};
        craux1 = reshape(craux1, ry(i)*n(i)*ry(i+1), 1);
        cry1 = cry1+auxfun(craux1);
    end;
    
    if (i>1)        
        cry1 = reshape(cry1, ry(i), n(i)*ry(i+1));
        [v,ilocr{i}] = max(abs(cry1.'));
        ilocr{i} = ilocr{i}(:);
        mid_inv{i,1} = 1/cry1(:,ilocr{i});
        mid_inv{i,2} = 1;
        Jyr{i} = indexmerge((1:n(i))', Jyr{i+1});
        Jyr{i} = Jyr{i}(ilocr{i},:);
        
        if (~isempty(aux))
            phiauxr{i} = reshape(aux{i}, raux(i)*n(i), raux(i+1));
            phiauxr{i} = phiauxr{i}*phiauxr{i+1};
            phiauxr{i} = reshape(phiauxr{i}, raux(i), n(i)*ry(i+1));
            phiauxr{i} = phiauxr{i}(:, ilocr{i});
        end;
    end;
    
    y{i} = reshape(cry1, ry(i), n(i), ry(i+1));
end;
for i=1:d
    J = indexmerge(Jyl{i}, (1:n(i))', Jyr{i+1});
    evalcnt = evalcnt + size(J,1);
    cry1 = autovecfun(fun, J, vec);
    if (~isempty(aux))
        craux1 = reshape(aux{i}, raux(i), n(i)*raux(i+1));
        craux1 = phiauxl{i}*craux1;
        craux1 = reshape(craux1, ry(i)*n(i), raux(i+1));
        craux1 = craux1*phiauxr{i+1};
        craux1 = reshape(craux1, ry(i)*n(i)*ry(i+1), 1);
        cry1 = cry1+auxfun(craux1);
    end;
    
    if (i<d)
        cry1 = reshape(cry1, ry(i)*n(i), ry(i+1));
        [v,ilocl{i+1}] = max(abs(cry1));
        mid_inv{i+1,1} = 1/cry1(ilocl{i+1});
        mid_inv{i+1,2} = 1;        
        Jyl{i+1} = indexmerge(Jyl{i}, (1:n(i))');
        Jyl{i+1} = Jyl{i+1}(ilocl{i+1},:);
        if (~isempty(aux))
            phiauxl{i+1} = reshape(aux{i}, raux(i), n(i)*raux(i+1));
            phiauxl{i+1} = phiauxl{i}*phiauxl{i+1};
            phiauxl{i+1} = reshape(phiauxl{i+1}, ry(i)*n(i), raux(i+1));
            phiauxl{i+1} = phiauxl{i+1}(ilocl{i+1}, :);
        end;
    end;
    
    y{i} = reshape(cry1, ry(i), n(i), ry(i+1));
end;

last_sweep = false;
maxy = 0;
max_dx = 0;
swp = 1;
dir = 1;
i = 1;
while (swp<=nswp)
    % Generate a random test set for new indices
    % Prepare candidate index sets--the ones with the current crosses excluded 
    cind1 = (1:ry(i)*n(i))';
    cind2 = (1:n(i+1)*ry(i+2))';
    cind1(ilocl{i+1})=[];
    cind2(ilocr{i+1})=[];
    
    if (strcmp(locsearch, 'als'))
        testsz = min(numel(cind1), numel(cind2));
    else
        %%% % rn Lottery
        % Now draw random entries
        tind = rand(min(numel(cind1), numel(cind2)), 1)*numel(cind1)*numel(cind2);
        tind = round(tind);
        tind(tind>numel(cind1)*numel(cind2))=numel(cind1)*numel(cind2);
        tind(tind<1)=1;
        tind = unique(tind);
        testsz = size(tind,1);
    end;
    
    if (~isempty(aux))
        craux1 = reshape(aux{i}, raux(i), n(i)*raux(i+1));
        craux1 = phiauxl{i}*craux1;
        craux1 = reshape(craux1, ry(i)*n(i), raux(i+1));
        craux2 = reshape(aux{i+1}, raux(i+1)*n(i+1), raux(i+2));
        craux2 = craux2*phiauxr{i+2};
        craux2 = reshape(craux2, raux(i+1), n(i+1)*ry(i+2));
    end;
    
    % Check that we are not in the full rank case
    if (testsz>0)
        if (strcmp(locsearch, 'als'))
            % 2D ALS cross
            % Evaluate y at tind
            cry1 = reshape(y{i}, ry(i)*n(i), ry(i+1));
            cry2 = reshape(y{i+1}, ry(i+1), n(i+1)*ry(i+2));
            ys1 = cry1*mid_inv{i+1,1};
            ys1 = ys1(cind1,:);
            ys2 = mid_inv{i+1,2}*cry2;
            ys2 = ys2(:,cind2);
            % Full indices
            J1 = indexmerge(Jyl{i}, (1:n(i))');
            J1c = J1(cind1,:);
            J2 = indexmerge((1:n(i+1))', Jyr{i+2});
            J2c = J2(cind2,:);
            rz = 2;
            cre2 = randn(numel(cind2), rz);
            [cre2,rv]=qr(cre2,0);
            indr = maxvol2(cre2);
            Jr = J2c(indr,:);
            Ye2 = ys2(:,indr);
            
            J = indexmerge(J1c,Jr);
            evalcnt = evalcnt + size(J,1);
            cre1 = autovecfun(fun, J, vec);
            if (~isempty(aux))
                craux = craux1(cind1,:)*craux2(:,cind2(indr));
                craux = craux(:);
                craux = auxfun(craux);
                cre1 = cre1+craux;
            end;
            maxy = max(maxy, max(abs(cre1)));
            cre1 = reshape(cre1, numel(cind1), rz);
            cre1 = cre1-ys1*Ye2;
            [zmax1,imaxnew]=max(abs(cre1(:)));
            imaxnew = tt_ind2sub([numel(cind1), rz], imaxnew);
            imax1 = imaxnew(1);
            [cre1,rv]=qr(cre1,0);
            indl = maxvol2(cre1);
            Jl = J1c(indl,:);
            Ye1 = ys1(indl,:);
            
            J = indexmerge(Jl,J2c);
            evalcnt = evalcnt + size(J,1);
            cre2 = autovecfun(fun, J, vec);
            if (~isempty(aux))
                craux = craux1(cind1(indl),:)*craux2(:,cind2);
                craux = craux(:);
                craux = auxfun(craux);
                cre2 = cre2+craux;
            end;
            maxy = max(maxy, max(abs(cre2)));
            cre2 = reshape(cre2, rz, numel(cind2));
            cre2 = cre2-Ye1*ys2;
            [zmax2,imaxnew]=max(abs(cre2(:)));
            imaxnew = tt_ind2sub([rz, numel(cind2)], imaxnew);
            imax2 = imaxnew(2);
            imax1 = cind1(imax1);
            imax2 = cind2(imax2);
            emax = max(zmax1,zmax2);
        else
            %%% % rn Lottery
            tind = tt_ind2sub([numel(cind1), numel(cind2)], tind);
            tind = [cind1(tind(:,1)), cind2(tind(:,2))];
            % Full indices
            J1 = indexmerge(Jyl{i}, (1:n(i))');
            J1c = J1(tind(:,1),:);
            J2 = indexmerge((1:n(i+1))', Jyr{i+2});
            J2c = J2(tind(:,2),:);
            J = [J1c,J2c];
            % Evaluate the function
            evalcnt = evalcnt + testsz;
            crt = autovecfun(fun, J, vec);
            if (~isempty(aux))
                craux = zeros(testsz, 1);
                for j=1:testsz
                    craux(j) = craux1(tind(j,1),:)*craux2(:,tind(j,2));
                end;
                crt = crt+auxfun(craux);
            end;
            maxy = max(maxy, max(abs(crt)));
            % Evaluate y at tind
            cry1 = reshape(y{i}, ry(i)*n(i), ry(i+1));
            cry2 = reshape(y{i+1}, ry(i+1), n(i+1)*ry(i+2));
            % Subtract the current approx.
            cry = zeros(testsz, 1);
            % Apply the inverse interp at the middle
            cre1 = cry1*mid_inv{i+1,1};
            cre2 = mid_inv{i+1,2}*cry2;
            for j=1:testsz
                cry(j) = cre1(tind(j,1),:)*cre2(:,tind(j,2));
            end;
            cre = crt-cry;
            % Take the max-err index to enrich
            [emax,imax2] = max(abs(cre));
            
            % All again: at (:,tind(imax,2)), run the error maximization
            J1c = J1(cind1,:);
            J = indexmerge(J1c, J2c(imax2,:));
            evalcnt = evalcnt + size(J,1);
            crt = autovecfun(fun, J, vec);
            if (~isempty(aux))
                craux = craux1(cind1,:)*craux2(:,tind(imax2,2));
                crt = crt+auxfun(craux);
            end;
            maxy = max(maxy, max(abs(crt)));
            % Subtract the current approx.
            cry = cre1(cind1,:)*cre2(:,tind(imax2,2));
            cre = crt-cry;
            % Take the max-err index to enrich
            [emax,imax1] = max(abs(cre));
            imax1 = cind1(imax1);
            imax2 = tind(imax2,2);
            % imax1 samples from ry(i)*n(i)
            % imax2 samples from n(i+1)*ry(i+2)
        end;
        
        dx = emax/maxy;
        max_dx = max(max_dx,dx);
        if (verb>1)
            fprintf('=greedy_cross= i=%d, swp=%d, testsz=%d, emax=%3.3e, dx=%3.3e, cond=[%3.3e,%3.3e]\n', i, swp, testsz, emax, dx, cond(mid_inv{i+1,1}), cond(mid_inv{i+1,2}));
        end;
        
        if (dx>tol)
            % Now evaluate new factors
            J1m = J1(imax1,:);
            J2m = J2(imax2,:);
            % Generate a full index = [Jl, (1:n)', Jr]
            % Jyl to lowest index           n(i) to the middle one,                      Jyr to the senior index
            Jl = indexmerge(J1, J2m);
            Jr = indexmerge(J1m, J2);
            
            evalcnt = evalcnt + size(Jl,1);
            cre1 = autovecfun(fun, Jl, vec);
            if (~isempty(aux))
                craux = craux1*craux2(:,imax2);
                cre1 = cre1+auxfun(craux);
            end;
            
            evalcnt = evalcnt + size(Jr,1);
            cre2 = autovecfun(fun, Jr, vec);
            if (~isempty(aux))
                craux = craux1(imax1,:)*craux2;
                craux = craux(:);
                cre2 = cre2+auxfun(craux);
            end;            
            
            cre1 = reshape(cre1, ry(i)*n(i), 1);
            cre2 = reshape(cre2, 1, n(i+1)*ry(i+2));
            
            % Expand blocks
            y{i} = [cry1, cre1];
            y{i+1} = [cry2; cre2];
            ry(i+1)=ry(i+1)+1;
            % ilocl{i+2} is smashed now, since ry(i+1) changed. Recover...
            if (i<d-1)
                ii = floor((ilocl{i+2}-1)/(ry(i+1)-1));
                ilocl{i+2} = ilocl{i+2}+ii;
            end;
            % Expand indices
            Jyl{i+1} = [Jyl{i+1}; J1m];
            Jyr{i+1} = [Jyr{i+1}; J2m];
            ilocl{i+1} = [ilocl{i+1}; imax1];            
            ilocr{i+1} = [ilocr{i+1}; imax2];
            % We know that the inverse LU may be written analytically
            uold = mid_inv{i+1,1};
            lold = mid_inv{i+1,2};
            % Expanding vectors
            erow = cry1(imax1, :);
            ecol = cre1(ilocl{i+1}(1:ry(i+1)-1), 1);
            % ..and the scalar
            eel = cre1(imax1,1);
            ecol = lold*ecol;
            erow = erow*uold;
            eel = eel - erow*ecol; % alpha-dA^{-1}c
            ecol = uold*ecol; % A^{-1} c
            erow = erow*lold; % d A^{-1}
            % New inv(U)
            mid_inv{i+1,1} = zeros(ry(i+1), ry(i+1));
            mid_inv{i+1,1}(1:ry(i+1)-1, 1:ry(i+1)-1)=uold;
            mid_inv{i+1,1}(1:ry(i+1)-1, ry(i+1)) = -ecol/eel;
            mid_inv{i+1,1}(ry(i+1), ry(i+1))=1/eel;
            % New inv(L)
            mid_inv{i+1,2} = zeros(ry(i+1), ry(i+1));
            mid_inv{i+1,2}(1:ry(i+1)-1, 1:ry(i+1)-1)=lold;
            mid_inv{i+1,2}(ry(i+1), 1:ry(i+1)-1) = -erow;
            mid_inv{i+1,2}(ry(i+1), ry(i+1))=1;
            
            y{i} = reshape(y{i}, ry(i), n(i), ry(i+1));
            y{i+1} = reshape(y{i+1}, ry(i+1), n(i+1), ry(i+2));                
        end;
    end;         
    
    if (~isempty(aux))
        phiauxl{i+1} = craux1(ilocl{i+1}, :);
        phiauxr{i+1} = craux2(:, ilocr{i+1});
    end;
    
    i = i+dir;
    
    % Check the convergence, restart
    if ((i==d)||(i==0))
        if (verb>0)
            fprintf('=greedy_cross= swp=%d, max_dx=%3.3e, max_rank=%d, cum#evals=%d\n', swp, max_dx, max(ry), evalcnt);
        end;
        
        if (dir>0)&&(last_sweep)
            break;
        end;
        
        if (max_dx<tol_exit)
            last_sweep = true;
        end;
        
        max_dx = 0;
        i=1;
        swp = swp+1;
    end;
end;

% Merge mid_inv: y:=inv(L)*y*inv(U)
for i=1:d
    y{i} = reshape(y{i}, ry(i), n(i)*ry(i+1));
    y{i} = mid_inv{i,2}*y{i};
    y{i} = reshape(y{i}, ry(i)*n(i), ry(i+1));
    y{i} = y{i}*mid_inv{i+1,1};
    y{i} = reshape(y{i}, ry(i), n(i), ry(i+1));
end;
y = cell2core(tt_tensor,y);
end


function [y]=autovecfun(fun, J, vec)
if (vec)
    y = fun(J); % User can cast M x d -> M x 1
else
    % They can't -- run a loop manually
    sz = size(J,1);
    y = zeros(sz, 1);
    for j=1:sz
        y(j) = fun(J(j,:));
    end;
end;
end


function [J]=indexmerge(varargin)
% Merges two or three indices in the little-endian manner
sz1 = max(size(varargin{1},1),1);
sz2 = max(size(varargin{2},1),1);
sz3 = 1;
if (nargin>2) % Currently allows only 3
    sz3 = max(size(varargin{3}, 1), 1);
end;
% J1 goes to the fastest index, just copy it
J1 = repmat(varargin{1}, sz2*sz3, 1);
% J2 goes to the middle
J2 = reshape(varargin{2}, 1, []);
J2 = repmat(J2, sz1, 1); % now sz1 ones will be the fastest
J2 = reshape(J2, sz1*sz2, []);
J2 = repmat(J2, sz3, 1);
J = [J1,J2];
if (nargin>2)
    % J3 goes to the slowest
    J3 = reshape(varargin{3}, 1, []);
    J3 = repmat(J3, sz1*sz2, 1); % now sz1 ones will be the fastest
    J3 = reshape(J3, sz1*sz2*sz3, []);
    J = [J,J3];
end;
end


