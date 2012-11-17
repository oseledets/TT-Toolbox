function [x] = tt_meshgrid(varargin)
%Analogue of the meshgrid function for the TT-format
%   X = TT_MESHGRID(A,B,C,...) Computes the meshgrid based on "1d"
%       representations 
%   X = TT_MESHGRID(A) Computes the meshgrid based on the cell array A of
%       the representations
%   X = TT_MESHGRID(T,D) Computes the d-dimensional meshgrid, using T as a
%       one-dimensional grid
if ( numel(varargin) == 2 )
    if ( ismatrix(varargin{2}) )
        d = varargin{2};
        t = varargin{1};
        z = cell(d,1);
        for i = 1:d
            z{i} = t;
        end
    else
        z = varargin{1};
    end
else
    z = varargin{1}; 
end
if ~iscell(z)
    z = varargin;
end
d = numel(z);
x = cell(d,1);
e = cell(d,1);
for i = 1:d
    e{i} = tt_ones(size(z{i}));
end
for i = 1:d
    v = z{i};
    for j = i-1:-1:1
        v = kron(e{j},v);
    end
    for j = i+1:d
        v = kron(v,e{j});
    end
    x{i} = v;
end
end