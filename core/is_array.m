function bool = is_array(x)
    % Check if x is regular dense array.
    bool = isa(x, 'double') || isa(x, 'single') || isa(x, 'gpuArray');
end
