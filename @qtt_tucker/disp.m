function disp(tt,name)
%DISP Command window display of a qtt_tucker.
%

if ~exist('name','var')
    name = 'ans';
end
fprintf('%s is a QTT-tucker tensor: \n',name);
fprintf('Tucker ranks: \n');
sz=size(tt.core);
rk=rank(tt.core);
for i=1:tt.dphys-1
   fprintf('%d-',sz(i)); 
end
fprintf('%d\n',sz(tt.dphys));
fprintf('TT-ranks for the core: \n');
for i=1:tt.dphys
   fprintf('%d-',rk(i)); 
end
fprintf('%d\n',rk(tt.dphys+1));
tk=tt.tuck;
fprintf('TT-ranks for the factors: \n');
for i=1:tt.dphys
   di=ndims(tk{i});
   ri=rank(tk{i});
   fprintf('Factor %d: ',i);
   for j=1:di-1
      fprintf('%d-',ri(j));
   end
   fprintf('%d\n',ri(di));
end
% %fprintf('r(1)=%d \n', r(1));
% for i=1:d
%    fprintf('r(%d)=%d \t n(%d)=%d \n',i,r(i),i,n(i));
% end
% fprintf('r(%d)=%d \n',d+1,r(d+1));

