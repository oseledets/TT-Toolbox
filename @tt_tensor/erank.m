function [er]=erank(tt)
%[ER]=ERANK(TT)
%Effective rank of the TT-tensor
n=tt.n; d=tt.d; r=tt.r;
sz=dot(n.*r(1:d),r(2:d+1));
er=sz./sum(n); er=sqrt(er);
return
end