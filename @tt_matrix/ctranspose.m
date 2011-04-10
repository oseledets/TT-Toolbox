function [tt1]=ctranspose(tt)
%[TT1]=CTRANSPOSE(TT)
%Compute complex conjugate transpose of a TT-matrix

%Zaglushka
tt1=tt;
tt1.tt=tt_tensor(tt_mat_to_vec(tt_transp(tt_vec_to_mat(core(tt.tt),tt.n,tt.m))));
%keyboard