function disp(tt,name)
%DISP Command window display of a tt_array.
%

if ~exist('name','var')
    name = 'ans';
end
tt1=tt.tt;
r=tt1.r; ns=tt.ns;
d=tt1.d;
fprintf('%s is a %dx%d-dimensional TT-array, ranks and mode sizes: \n',name,d,size(ns,2));
%fprintf('r(1)=%d \n', r(1));
for i=1:d
    curstr = sprintf('r(%d)=%d \t ', i, r(i));
    for j=1:size(ns,2)
        curstr = [curstr, sprintf('%d  ', ns(i,j))];
    end;
    fprintf([curstr, '\n']);
%    fprintf('r(%d)=%d \t n(%d)=%d m(%d) = %d \n',i,r(i),i,n(i),i,m(i));
end
fprintf('r(%d)=%d \n',d+1,r(d+1));

