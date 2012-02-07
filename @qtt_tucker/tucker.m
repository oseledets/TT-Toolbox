function [factors,cr]=tucker(tt,eps)
%[factors,cr]=tucker(tt,eps)
%Get Tucker factors and Tucker core in the TT-format
%for a given TT-tensor

%Zaglushka
[factors,cr]=tt_tuck(tt,eps);
return
end