function disp(tt,name,varargin)
%Command window display of a QTT_TUCKER
%
disp_mode_sizes=false;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'modes'
            disp_mode_sizes=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

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
fprintf('TT-ranks/mode sizes for the factors: \n');
for i=1:tt.dphys
   di=ndims(tk{i});
   ri=rank(tk{i});
   ni=size(tk{i});
   fprintf('Factor %d: ',i);
   for j=1:di
      fprintf('%d-',ri(j));
   end   
   fprintf('%d\n',ri(di+1));
   if ( disp_mode_sizes ) 
     fprintf('Modes   : ');
     for j=1:di-1
       fprintf('%d-',ni(j));
     end
     fprintf('%d\n',ni(di));
   end
end
% %fprintf('r(1)=%d \n', r(1));
% for i=1:d
%    fprintf('r(%d)=%d \t n(%d)=%d \n',i,r(i),i,n(i));
% end
% fprintf('r(%d)=%d \n',d+1,r(d+1));

