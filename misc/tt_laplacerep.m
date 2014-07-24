% Creates the Laplace-like operator S x M x M + M x S x M + M x M x S in
% the case when S and M are given as tt_matrices themselves (e.g. QTT or
% spin operators). S and M can be given as d x 1 cell arrays, or single
% tt_matrices. In the latter case, the third input d must be specified.

function [A] = tt_laplacerep(S,M,d)

% Grumble inputs
if (nargin<3)
    d = [];
end;
ds = [];
if (isa(S,'cell'))
    ds = numel(S);
end;
dm = [];
if (isa(M,'cell'))
    dm = numel(M);
end;
d = min([ds,dm,d]);
if (isempty(d))
    error('No replication dimension is given');
end;

% Check if no replication is needed
if (d==1)
    A = S;
    if (isa(A,'cell'))
        A = A{1};
    end;
    return;
end

nm = [];
if (isa(S, 'tt_matrix'))
    S1 = S;
    S = cell(d,1);
    for i=1:d
        S{i} = S1;
        nm = [nm; [S1.n, S1.m]];
    end;
else
    for i=1:d
        nm = [nm; [S{i}.n, S{i}.m]];
    end;
end;
if (isa(M, 'tt_matrix'))
    M1 = M;
    M = cell(d,1);
    for i=1:d
        M{i} = M1;
    end;
end;

A = cell(size(nm,1), 1);

% First "block". Consisting of S{1}.d blocks in fact..
rs = S{1}.r;
rm = M{1}.r;
% First first block
A{1} = zeros(1, nm(1,1), nm(1,2), rs(2)+rm(2));
A{1}(1,:,:,1:rs(2)) = S{1}{1};
A{1}(1,:,:,rs(2)+1:rs(2)+rm(2)) = M{1}{1};
% Rest first blocks
for i=2:S{1}.d
    A{i} = zeros(rs(i)+rm(i), nm(1,1), nm(1,2), rs(i+1)+rm(i+1));
    A{i}(1:rs(i),:,:,1:rs(i+1)) = S{1}{i};
    A{i}(rs(i)+1:rs(i)+rm(i),:,:,rs(i+1)+1:rs(i+1)+rm(i+1)) = M{1}{i};
end;

% The global dimension index
i_glob = S{1}.d+1;

% Middle blocks
for j=2:d-1
    L = S{j}.d;
    if (L<2)
        error('dim(S{%d})<2 is not supported', j);
    end;
    rs = S{j}.r;
    rm = M{j}.r;
    % First middle block
    A{i_glob} = zeros(2, nm(i_glob,1), nm(i_glob,2), rm(2)+rs(2)+rm(2));
    A{i_glob}(1,:,:,1:rm(2)) = M{j}{1};
    A{i_glob}(2,:,:,rm(2)+1:rm(2)+rs(2)) = S{j}{1};
    A{i_glob}(2,:,:,rm(2)+rs(2)+1:rm(2)+rs(2)+rm(2)) = M{j}{1};
    i_glob = i_glob+1;
    % Middle Middle blocks
    for i=2:L-1
        A{i_glob} = zeros(rm(i)+rs(i)+rm(i), nm(i_glob,1), nm(i_glob,2), rm(i+1)+rs(i+1)+rm(i+1));
        A{i_glob}(1:rm(i),:,:,1:rm(i+1)) = M{j}{i};
        A{i_glob}(rm(i)+1:rm(i)+rs(i),:,:,rm(i+1)+1:rm(i+1)+rs(i+1)) = S{j}{i};
        A{i_glob}(rm(i)+rs(i)+1:rm(i)+rs(i)+rm(i),:,:,rm(i+1)+rs(i+1)+1:rm(i+1)+rs(i+1)+rm(i+1)) = M{j}{i};
        i_glob = i_glob+1;
    end;
    % Last middle block
    A{i_glob} = zeros(rm(L)+rs(L)+rm(L), nm(i_glob,1), nm(i_glob,2), 2);
    A{i_glob}(1:rm(L),:,:,1) = M{j}{L};
    A{i_glob}(rm(L)+1:rm(L)+rs(L),:,:,1) = S{j}{L};
    A{i_glob}(rm(L)+rs(L)+1:rm(L)+rs(L)+rm(L),:,:,2) = M{j}{L};
    i_glob = i_glob+1;
end;

% middle Last blocks
L = S{d}.d;
rs = S{d}.r;
rm = M{d}.r;
for i=1:L-1
    A{i_glob} = zeros(rm(i)+rs(i), nm(i_glob,1), nm(i_glob,2), rm(i+1)+rs(i+1));
    A{i_glob}(1:rm(i),:,:,1:rm(i+1)) = M{d}{i};
    A{i_glob}(rm(i)+1:rm(i)+rs(i),:,:,rm(i+1)+1:rm(i+1)+rs(i+1)) = S{d}{i};
    i_glob = i_glob+1;
end;
% Last Last block
A{i_glob} = zeros(rm(L)+rs(L), nm(i_glob,1), nm(i_glob,2), 1);
A{i_glob}(1:rm(L),:,:,1) = M{d}{L};
A{i_glob}(rm(L)+1:rm(L)+rs(L),:,:,1) = S{d}{L};

A = cell2core(tt_matrix,A);

end