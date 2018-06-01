function [X]=tt_meshgrid_vert(varargin)
% Analogue of the meshgrid function for the TT-format,
% !!!!!!!!!!!! concatenating tensor trains with ones with tkron !!!!!!!!!!!
% !!!!!!!!!!!! contrarily to tt_meshgrid (which uses kron2)     !!!!!!!!!!!
%   X = TT_MESHGRID_VERT(A,B,C,...) Computes the meshgrid based on "1d"
%       representations 
%   X = TT_MESHGRID_VERT(A) Computes the meshgrid based on the cell array A of
%       the representations
%   X = TT_MESHGRID_VERT(T,D) Computes the d-dimensional meshgrid, using T as a
%       one-dimensional grid
%   X = TT_MESHGRID_VERT(..., 'id', M) Expands with M instead of all-ones

vars = varargin;
M = [];
for i=1:numel(vars)
    if (isa(vars{i}, 'char'))&&(strcmp(vars{i},'id'))
        M = vars{i+1};
        vars(i:i+1) = [];
        break;
    end
end

if (numel(vars)==1)&&(isa(vars{1}, 'cell'))
    % We have a cell array of 1d points
    X = vars{1};
    d = numel(X);
    X = reshape(X, 1, d);
elseif (numel(vars)==2)&&(isa(vars{1}, 'tt_tensor')||isa(vars{1}, 'tt_matrix'))&&(isscalar(vars{2}))
    % First argument is a single x, second argument is the dimension
    d = vars{2};
    X = cell(1, d);
    for i=1:d
        X{i} = vars{1}; % make the same format as in the first case
    end
else
    % A set of independent x
    d = numel(vars);
    X = cell(1,d);
    for i=1:d
        if (isa(vars{i}, 'tt_tensor')||isa(vars{i}, 'tt_matrix'))
            X{i} = vars{i};
        else
            error('wrong input %d to tt_meshgrid_vert', i);
        end
    end
end

% Now X contains 1D tt_tensors. Expand them by ones

% Concat all n to the common storage
n = zeros(0,1);
% positions in n where each component starts
pos = ones(d+1,1);
for i=1:d
    n = [n; X{i}.n];
    pos(i+1) = pos(i)+X{i}.d;
end
% Expand
for i=1:d
    if (i>1)
        if (isempty(M))
            expand = tt_ones(n(1:pos(i)-1));
        else
            expand = M;
            if (i>2)
                expand = mtkron(repmat({M},1,i-1));
            end
        end
        X{i} = tkron(expand, X{i});
    end
    if (i<d)
        if (isempty(M))
            expand = tt_ones(n(pos(i+1):pos(d+1)-1));
        else
            expand = M;
            if (i<d-1)
                expand = mtkron(repmat({M},1,d-i));
            end
        end        
        X{i} = tkron(X{i}, expand);
    end
end

end
