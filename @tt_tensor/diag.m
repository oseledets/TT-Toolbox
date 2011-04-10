function [tm]=diag(tt)
%[TM]=DIAG(TT)
%Constructs diagonal matrix from TT-tensor
%Zaglushka
tm=tt_matrix(tt_vec_to_diag(core(tt)));
return
end