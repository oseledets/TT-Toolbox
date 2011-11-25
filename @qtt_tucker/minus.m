function [a]=minus(b,c)
%[A]=MINUS(B,C)
%Subtract two QTT_Tucker tensors
a=b+(-1.0)*c;
return
end