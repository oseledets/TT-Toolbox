function x = move2gpu(x, toSingle)
    % Move TT-tensor to the GPU.
    % toSingle is boolean variable. If true, store everything on the GPU as
    % single.
    error('I''m pretty sure you don''t need this function. cnn_train calls vl_simplenn_move, which moves all layers (including the TT ones) to GPU when appropriate.')
    if nargin <= 1
        toSingle = false;
    end
    if toSingle
        x.core = single(x.core);
    end
    x.core = gpuArray(x.core);
end
