function [mem] = tt_save(filename, tt, precision)
% function [mem] = tt_save(filename, tt, [precision])
% Saves a tt_tensor tt to a filename in SDV =) format 
% precision: 0 - double, 1 - single
% returns the bytes written to a file,
% -1 in case of errors

if (nargin<3)||(isempty(precision))
    precision = 0; % double
end;

mem = tt_write(filename, tt.d, tt.r, tt.n, tt.core, precision);

end
