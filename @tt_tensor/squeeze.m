function [tt]=squeeze(tt)
%[TT]=SQUEEZE(TT)
%Removes singleton dimensions from the TT-tensor

%Zaglushka
tt=tt_tensor(tt_squeeze(core(tt)));
return
end