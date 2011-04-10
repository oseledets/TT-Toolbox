% FUNCTION [TTM]=TT_DIAG(TT_A)
% Returns diagonal matrix in QTTM format with the values of vector 
% tt_a in QTT format along the diagonal

function [ttm] = tt_diag(tt_a)

d = size(tt_a, 1);
ttm = cell(d,1);
for i=1:d
    ttm{i}=zeros(2,2, size(tt_a{i},2), size(tt_a{i},3));
    ttm{i}(1,1,:,:) = reshape(tt_a{i}(1,:,:), 1,1,size(tt_a{i},2), size(tt_a{i},3));
    ttm{i}(2,2,:,:) = reshape(tt_a{i}(2,:,:), 1,1,size(tt_a{i},2), size(tt_a{i},3));
end;

end