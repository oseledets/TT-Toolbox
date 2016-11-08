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

if (nargin==1)&&(isa(varargin{1}, 'cell'))
    % We have a cell array of 1d points
    X = varargin{1};
    d = numel(X);
    X = reshape(X, 1, d);
elseif (nargin==2)&&(isa(varargin{1}, 'tt_tensor'))&&(isscalar(varargin{2}))
    % First argument is a single x, second argument is the dimension
    d = varargin{2};
    X = cell(1, d);
    for i=1:d
        X{i} = varargin{1}; % make the same format as in the first case
    end;
else
    % A set of independent x
    d = nargin;
    X = cell(1,d);
    for i=1:d
        if (isa(varargin{i}, 'tt_tensor'))
            X{i} = varargin{i};
        else
            error('wrong input %d to tt_meshgrid_vert', i);
        end;
    end;
end;

% Now X contains 1D tt_tensors. Expand them by ones

% Concat all n to the common storage
n = zeros(0,1);
% positions in n where each component starts
pos = ones(d+1,1);
for i=1:d
    n = [n; X{i}.n];
    pos(i+1) = pos(i)+X{i}.d;
end;
% Expand
for i=1:d
    if (i>1)
        X{i} = tkron(tt_ones(n(1:pos(i)-1)), X{i});
    end;
    if (i<d)
        X{i} = tkron(X{i}, tt_ones(n(pos(i+1):pos(d+1)-1)));
    end;
end;

end