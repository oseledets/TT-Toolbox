function enlarged_ttm = add_non_essential_dims(small_ttm, new_n, old_vars_pos)
    % enlarged_ttm = add_non_essential_dims(small_ttm, new_n, old_vars_pos)
    % creates enlarged version of tt_matrix, expanding it by identities.
    % If you need to expand a matrix by ones, use add_non_essential_dims(small_ttm.tt, ...).
    % Mode sizes of the enlarged_ttm are equal to new_n.
    % First dimension of the small_tt is old_vars_pos(1) dimension of the enlarged_ttm.
    % This procedure doesn't increase the maximal TT-rank.
    % 
    % Input:
    % 	small_ttm    -- P-dimensional matrix in the TT-format.
    % 	new_n        -- Q by 2 vector with enlarged_tt mode sizes (Q > P).
    % 	old_vars_pos -- P by 1 vector in ascending order. Consists of positions
    % 							of the small_ttm dimensions in the enlarged_ttm tensor.
    % 
    % Output:
    % 	enlarged_ttm  -- Q-dimensional matrix in the TT format.
    % 			enlarged_ttm(i1,j1, ..., iQ,jQ)
    % 			i1 take values from {1, ..., new_n(1,1)}
    % 			j1 take values from {1, ..., new_n(1,2)}    
    % 			...
    % 			iQ take values from {1, ..., new_n(Q,1)}
    % 			jQ take values from {1, ..., new_n(Q,2)}    
    % 
    % Example:    
    %       small_ttm = tt_shift(2,3,2);
    %       enlarged_ttm = add_non_essential_dims(small_ttm, [2;3;2;3;2;3]*[1,1], [1;3;5]);
    %   enlarged_ttm is the shift in dimensions 1,3,5 and the identity in
    %   dimensions 2,4,6.
    

    if (~issorted(old_vars_pos))
        error('Dimensions (third argument) must be in ascending order.');
    end

    small_ttm = small_ttm.tt;
    
    if size(new_n,2)~=2 || ~isvector(old_vars_pos) || small_ttm.d ~= length(old_vars_pos) || size(new_n,1)<length(old_vars_pos)
    	error('Wrong usage, see help tt_matrix.add_non_essential_dims.');
    end



    d = size(new_n,1);
    enlarged_ttm = tt_tensor; % will convert to tt_matrix in the end
    enlarged_ttm.d = d;
    enlarged_ttm.r = zeros(d + 1, 1);
    enlarged_ttm.r(1) = 1;
    enlarged_ttm.n = prod(new_n, 2);
    enlarged_ttm.core = [];
    enlarged_ttm.ps = zeros(d + 1, 1);
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
            if new_n(dimension_i,1)*new_n(dimension_i,2) ~= small_ttm.n(curr_var_i)
                error('Dimensions must agree!');
            end
        end

        if before_first_var || after_last_var
            % Dimension_i is before old_vars_pos(1) or after old_vars_pos(end).
            if new_n(dimension_i,1)~=new_n(dimension_i,2)
                error('Row and column sizes in dimension %d must be the same', dimension_i);
            end;
            curr_core = eye(new_n(dimension_i,1));
            curr_rank = 1;
        elseif old_vars_pos(curr_var_i) == dimension_i
            % It's real dimension.
            curr_core = core(small_ttm, curr_var_i);
            curr_rank = small_ttm.r(curr_var_i);
        else
            % Non-essential dimension between two real ones.
            curr_rank = small_ttm.r(curr_var_i + 1);
            if new_n(dimension_i,1)~=new_n(dimension_i,2)
                error('Row and column sizes in dimension %d must be the same', dimension_i);
            end;            
            curr_core = core_eye(curr_rank, new_n(dimension_i,1));
        end
        enlarged_ttm.r(dimension_i) = curr_rank;
        enlarged_ttm.ps(dimension_i) = length(enlarged_ttm.core) + 1;
        enlarged_ttm.core = [enlarged_ttm.core; curr_core(:)];
    end
    enlarged_ttm.r(end) = 1;
    enlarged_ttm.ps(end) = length(enlarged_ttm.core) + 1;
    
    enlarged_ttm = tt_matrix(enlarged_ttm, new_n(:,1), new_n(:,2));
end

function core = core_eye(r, n)
    % Return 4-dimensional array filled with eye matrices,
    % which corresponds to non-essential dimension core.
    core = zeros(r, n, n, r);
    for i = 1:n
        core(:, i, i, :) = eye(r);
    end
end
