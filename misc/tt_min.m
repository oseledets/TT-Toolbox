function [res]=tt_min(tt)
% function [res]=tt_min(tt)
% Finds a (quasi)minimum of the TT-tensor tt, provided
% the maximum is positive

resmax = tt_max(tt);
ons = tt_ones(size(tt,1), tt_size(tt));

res = resmax - tt_max(tt_add(tt, tt_scal(ons,-resmax)));

end
