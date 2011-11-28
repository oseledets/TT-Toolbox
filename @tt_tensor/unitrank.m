function [tt]=unitrank(tt)
% function [tt]=unitrank(tt)
% Stuff the first and the last ranks to the first and the last modes, 
% making the boundary ranks to be equal to 1.
% The first mode represents now (\alpha_0 n_1), and the last (n_d \alpha_d)
% Probably would help to work with tts<->rakes

d = tt.d;
tt.n(1) = tt.r(1)*tt.n(1);
tt.r(1) = 1;
tt.n(d) = tt.n(d)*tt.r(d+1);
tt.r(d+1) = 1;

end