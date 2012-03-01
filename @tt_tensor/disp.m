function str=disp(tt,name)
%Command window display of a TT-tensor.
%
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

if ~exist('name','var')
    name = 'ans';
end
r=tt.r; n=tt.n; 
d=tt.d;
str=[];
str=[str,sprintf('%s is a %d-dimensional TT-tensor, ranks and mode sizes: \n',name,d)];
%fprintf('r(1)=%d \n', r(1));
for i=1:d
    str=[str,sprintf('r(%d)=%d \t n(%d)=%d \n',i,r(i),i,n(i))];
%   fprintf(' \t n(%d)=%d \n',i,n(i));
%   fprintf('r(%d)=%d \n',i+1,r(i+1));
end
 str=[str,sprintf('r(%d)=%d \n',d+1,r(d+1))];
 if ( nargout == 0 )
    fprintf(str);
 end
end

