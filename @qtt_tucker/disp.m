function disp(tt,name)
%DISP Command window display of a tt_tensor.
%

if ~exist('name','var')
    name = 'ans';
end
r=tt.r; n=tt.n; 
d=tt.d;
fprintf('%s is a %d-dimensional TT-tensor, ranks and mode sizes: \n',name,d);
%fprintf('r(1)=%d \n', r(1));
for i=1:d
   fprintf('r(%d)=%d \t n(%d)=%d \n',i,r(i),i,n(i));
end
fprintf('r(%d)=%d \n',d+1,r(d+1));

