function [mm]=mem(tt)
%[MM]=MEM(TT)
%Computes memory for TT-tensor (only for cores)
mm=dot(tt.n.*tt.r(1:tt.d),tt.r(2:tt.d+1));
return
end