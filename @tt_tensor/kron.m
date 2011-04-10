function [c] = kron(a,b)
%[C]=KRON(A,B)
%Kronecker product of two TT-tensors

%Zaglushka
c=tt_tensor(tt_kron(core(a),core(b)));
return
end
