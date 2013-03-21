% function [tts] = tt_sin_cos(d, alpha, type)
% Returns array with i-th element equals to sin(alpha*(i-1)) or cos(alpha*(i-1))
% in QTT format of size 2^d;
% type:
%     0 - sin
%     1 - cos

function [tts] = tt_sin_cos(d, alpha, type)

tts = cell(d,1);
tts{1} = zeros(2,2);
if (type==0)
    tts{1}(1:2,1) = [0 sin(alpha)]';
    tts{1}(1:2,2) = [1 cos(alpha)]';
else
    tts{1}(1:2,1) = [1 cos(alpha)]';
    tts{1}(1:2,2) = [0 -sin(alpha)]';    
end;

for q=2:d-1
    tts{q} = zeros(2,2,2);
    tts{q}(1:2,1,1)=[1 cos((2^(q-1))*alpha)]';
    tts{q}(1:2,2,1)=[0 sin((2^(q-1))*alpha)]';
    tts{q}(1:2,1,2)=[0 -sin((2^(q-1))*alpha)]';
    tts{q}(1:2,2,2)=[1 cos((2^(q-1))*alpha)]';
end;

tts{d} = zeros(2,2);
tts{d}(1:2,1)=[1 cos((2^(d-1))*alpha)]';
tts{d}(1:2,2)=[0 sin((2^(d-1))*alpha)]';

tts = tt_tensor(tts);

end