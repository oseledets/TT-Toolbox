function enlarged_tt = add_non_essential_dims(small_tt, new_n, old_vars_pos)
    % enlarged_tt = add_non_essential_dims(small_tt, new_n, old_vars_pos) creates enlarged version of TT-tensor.
    % Mode sizes of the enlarged_tt are equal to new_n.
    % First dimension of the small_tt is old_vars_pos(1) dimension of the enlarged_tt.
    % This procedure doesn't increase the maximal TT-rank.
    % 
    % Input:
    % 	small_tt     -- P-dimensional tensor in the TT-format.
    % 	new_n        -- Q by 1 vector with enlarged_tt mode sizes (Q > P).
    % 	old_vars_pos -- P by 1 vector in ascending order. Consists of positions
    % 							of the small_tt dimensions in the enlarged_tt tensor.
    % 
    % Output:
    % 	enlarged_tt  -- Q-dimensional tensor in the TT format.
    % 			enlarged_tt(i1, ..., iQ)
    % 			i1 take values from {1, ..., new_n(1)}
    % 			...
    % 			iQ take values from {1, ..., new_n(Q)}
    % 
    % Example:
    %       m = rand(3, 4);
    %       small_tt = tt_tensor(m); % small_tt(i2, i4)
    %       enlarged_tt = add_non_essential_dims(small_tt, [2, 3, 6, 4, 4], [2, 4]);
    %       % enlarged_tt(i1, i2, i3, i4, i5, i6) == small_tt(i2, i4) for any i1, i3, i5, i6
    % 
    

    if (~issorted(old_vars_pos))
        error('Dimensions (third argument) must be in ascending order.');
    end

    if ~isvector(new_n) | ~isvector(old_vars_pos) | small_tt.d ~= length(old_vars_pos)
    	error('Wrong usage, see help add_non_essential_dims.');
    end



    d = length(new_n);
    enlarged_tt = tt_tensor;
    enlarged_tt.d = d;
    enlarged_tt.r = zeros(d + 1, 1);
    enlarged_tt.r(1) = 1;
    enlarged_tt.n = new_n(:);
    enlarged_tt.core = [];
    enlarged_tt.ps = zeros(d + 1, 1);
    enlarged_tt.over = 0;
    curr_var_i = 0;
    vars_length = length(old_vars_pos);
    before_first_var = true;
    after_last_var = false;
    for dimension_i = 1:d
        if curr_var_i == vars_length
            after_last_var = true;
        end

        if curr_var_i < vars_length && old_vars_pos(curr_var_i + 1) == dimension_i
            curr_var_i = curr_var_i + 1;
            before_first_var = false;
            if new_n(dimension_i) ~= small_tt.n(curr_var_i)
                error('Dimensions must agree!');
            end
        end

        if before_first_var || after_last_var
            % Dimension_i is before old_vars_pos(1) or after old_vars_pos(end).
            curr_core = ones(new_n(dimension_i), 1);
            curr_rank = 1;
        elseif old_vars_pos(curr_var_i) == dimension_i
            % It's real dimension.
            curr_core = core(small_tt, curr_var_i);
            curr_rank = small_tt.r(curr_var_i);
        else
            % Non-essential dimension between two real ones.
            curr_rank = small_tt.r(curr_var_i + 1);
            curr_core = core_eye(curr_rank, new_n(dimension_i));
        end
        enlarged_tt.r(dimension_i) = curr_rank;
        enlarged_tt.ps(dimension_i) = length(enlarged_tt.core) + 1;
        enlarged_tt.core = [enlarged_tt.core; curr_core(:)];
    end
    enlarged_tt.r(end) = 1;
    enlarged_tt.ps(end) = length(enlarged_tt.core) + 1;
end

function core = core_eye(r, n)
    % Return 3-dimensional array filled with eye matrices,
    % which corresponds to non-essential dimension core.
    core = zeros(r, n, r);
    for i = 1:n
        core(:, i, :) = eye(r);
    end
end
