function x = move2gpu(x, toSingle)
    % Move TT-tensor to the GPU.
    % toSingle is boolean variable. If true, store everything on the GPU as
    % single.
    if nargin <= 1
        toSingle = false;
    end
    if toSingle
        x.core = single(x.core);
    end
    x.core = gpuArray(x.core);
end
