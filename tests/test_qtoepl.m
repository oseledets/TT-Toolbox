function test_qtoepl()
%TEST_QTOEPL Compare multilevel QTT Toeplitz matrices with direct indexing.
% From the toolbox root: setup; addpath('tests'); test_qtoepl
% All cases are small and deterministic; no compiled routines are needed.

% A hand-computable case fixes the lag sign independently of the helper.
x = tt_tensor(reshape([0; 10; 20; 30], [2 2]));
assert(norm(full(tt_qtoepl(x, 1))-[20 10; 30 20], 'fro') < 1e-12);

% Include previously working 1D/equal/palindromic cases and unequal axis
% depths in both orders. Complex nonsymmetric data expose transposition,
% conjugation and axis-order errors that an even real kernel could hide.
depths = {1, 2, [2 2], [1 2 1], [1 2], [2 1], [1 2 3], [3 2 1]};
for k = 1:numel(depths)
    d = depths{k};
    n = 2.^d;
    j = (1:prod(2*n)).';
    values = sin(j)+1i*cos(j.^2);
    x = tt_tensor(reshape(values, 2*ones(1, sum(d+1))), 1e-13);
    expected = direct_toeplitz(values, n);
    for orientation = 1:2
        if orientation == 1, argument = d; else, argument = d.'; end
        actual = full(tt_qtoepl(x, argument));
        relative_error = norm(actual-expected, 'fro')/norm(expected, 'fro');
        assert(relative_error < 1e-11, ...
            'Incorrect Toeplitz matrix for depths %s (orientation %d): %.3g.', ...
            mat2str(d), orientation, relative_error);
    end
end

% Exhaust the entire kernel basis on the smallest unequal grid (2-by-4).
% Each unit kernel must produce exactly one spatial-shift matrix. This
% checks the complete linear map, including all 11 unused padding entries.
d = [1 2]; n = 2.^d;
points = coordinates((0:prod(n)-1).', n);
kernel_points = coordinates((0:prod(2*n)-1).', 2*n);
zero_images = 0;
for k = 1:prod(2*n)
    values = zeros(prod(2*n), 1); values(k) = 1;
    x = tt_tensor(reshape(values, 2*ones(1, sum(d+1))), 1e-13);
    lag = kernel_points(k, :)-(n+1);
    expected = ones(prod(n));
    for axis = 1:numel(n)
        differences = bsxfun(@minus, points(:, axis), points(:, axis).');
        expected = expected .* (differences == lag(axis));
    end
    zero_images = zero_images + ~any(expected(:));
    assert(norm(full(tt_qtoepl(x, d))-expected, 'fro') < 1e-12, ...
        'Incorrect spatial shift for kernel basis entry %d.', k);
end
assert(zero_images == 11);
fprintf('PASS test_qtoepl: direct complex matrices, row/column depths and all 32 kernel basis vectors.\n');
end

function matrix = direct_toeplitz(values, n)
% Definition: A(i,j) = X(N(1)+1+i(1)-j(1), ..., N(D)+1+i(D)-j(D)).
% This oracle uses only ordinary array indexing, not QTT transformations.
points = coordinates((0:prod(n)-1).', n);
strides = [1 cumprod(2*n(1:end-1))];
matrix = complex(zeros(prod(n)));
for row = 1:prod(n)
    indices = bsxfun(@minus, points(row, :), points);
    indices = bsxfun(@plus, indices, n+1);
    linear = 1+(indices-1)*strides.';
    matrix(row, :) = values(linear).';
end
end

function points = coordinates(indices, n)
% Decode physical indices in MATLAB column-major order, first axis fastest.
points = zeros(numel(indices), numel(n));
for axis = 1:numel(n)
    points(:, axis) = mod(indices, n(axis))+1;
    indices = floor(indices/n(axis));
end
end
