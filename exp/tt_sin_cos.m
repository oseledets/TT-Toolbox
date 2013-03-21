% function [tts] = tt_sin_cos(d, alpha, type)
% Returns array with i-th element equals to sin(alpha*(i-1)) or cos(alpha*(i-1))
% in QTT format of size 2^d;
% type:
%     0 - sin
%     1 - cos

function [tts] = tt_sin_cos(d, alpha, type)

tts = cell(d,1);
tts{1} = zeros(1, 2, 2);
if (type==0)
    tts{1}(1, 1:2, 1) = [1 cos(alpha)]';
    tts{1}(1, 1:2, 2) = [0 sin(alpha)]';
else
    tts{1}(1, 1:2, 1) = [0 -sin(alpha)]';
    tts{1}(1, 1:2, 2) = [1 cos(alpha)]';    
end;

for q=2:d-1
    tts{q} = zeros(2,2,2);
    % core matrix for n_q = 0:
    tts{q}(1:2,1,1)=[1 0]';
    tts{q}(1:2,1,2)=[0 1]';
    % core matrix for n_q = 1:
    tts{q}(1:2,2,1)=[cos((2^(q-1))*alpha) sin((2^(q-1))*alpha)]';
    tts{q}(1:2,2,2)=[-sin((2^(q-1))*alpha) cos((2^(q-1))*alpha)]';
end;

tts{d} = zeros(2,2);
% core col for n_d = 0:
tts{d}(1:2,1)=[0 1]';
% core col for n_d = 1:
tts{d}(1:2,2)=[sin((2^(d-1))*alpha) cos((2^(d-1))*alpha)]';

tt1 = tt_tensor();
%tt1.d = d;
%tt1.ps = tt1.get_ps();
%tt1.n = 2 * ones(1, d);
%tt1.r = [1; 2 * ones(1, d - 1); 1];
tt1 = cell2core(tt1, tts);
tts = tt1;
%tts = tt_tensor(tts);

end
