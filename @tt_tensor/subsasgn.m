function t = subsasgn(t,s,b)
%Indexing of cores of TT-format: T{I}=B
%   T{I} = B Set i-th core to be equal to B
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
 switch s(1).type
     case '.'
         switch s(1).subs
             case 'n'
                 if (size(b,1)==1)
                     b=b.';
                 end;
                 if (numel(s)>1)
                     s = s(2:end);
                     t.n = subsasgn(t.n, s, b);
                 else
                     t.n = b;
                 end;
%                  if (numel(s)>1)&&(strcmp(s(2).type,'()'))
%                      if (size(s(2).subs{1},1)==1)
%                          s(2).subs{1}=s(2).subs{1}';
%                      end;
%                      t.n(s(2).subs{1})=b;
%                  else
%                      t.n=b;
%                  end;
             case 'r'
                 if (size(b,1)==1)
                     b=b.';
                 end;
                 if (numel(s)>1)&&(strcmp(s(2).type,'()'))
                     if (size(s(2).subs{1},1)==1)
                         s(2).subs{1}=s(2).subs{1}';
                     end;
                     t.r(s(2).subs{1})=b;
                 else
                     t.r=b;
                 end;
             case 'core'                 
                 if (size(b,1)==1)
                     b=b.';
                 end
                 % I don't want to dispatch all the rest subs to the core,
                 % because MATLAB will copy it. Fortunately, only () is
                 % reasonable here.
                 if (numel(s)>1)&&(strcmp(s(2).type,'()'))
                     if (size(s(2).subs{1},1)==1)
                         s(2).subs{1}=s(2).subs{1}';
                     end;
                     t.core(s(2).subs{1})=b;
                 else
                     t.core=b;
                 end;
             case 'd'
                 t.d=b;
             case 'ps'
                 if (size(b,1)==1)
                     b=b.';
                 end
                 if (numel(s)>1)&&(strcmp(s(2).type,'()'))
                     if (size(s(2).subs{1},1)==1)
                         s(2).subs{1}=s(2).subs{1}';
                     end;
                     t.ps(s(2).subs{1})=b;
                 else
                     t.ps=b;
                 end
             otherwise
                 error(['Cannot change field ', s.subs, ' directly.']);
         end
     case '{}'
         % Insert one core (r1-n-r2!!!!!! not n-r1-r2!) into the TT2
         r1 = size(b,1);
         n1 = size(b,2);
         r2 = size(b,3);
         i = s.subs{1};
%          if (i==1)
%              r2 = r1;
%              r1 = 1;
%          end;
         n = t.n;
         d = t.d;
         r = t.r;
         cr = t.core;
         ps = t.ps;
         % Update the positions of the untouched cores
         psnew = ps;
         psnew(i+1:end) = ps(i+1:end)-(ps(i+1)-ps(i))+n1*r1*r2;
         cr(psnew(i+1):psnew(d+1)-1) = cr(ps(i+1):ps(d+1)-1);
%          if (i~=1)
%              b = permute(b, [2, 1, 3]); % We have r1,n1,r2 in the new format
%          end;
         % Update the current core
         cr(psnew(i):psnew(i+1)-1)=b(:);
         cr = cr(1:psnew(d+1)-1);
         % Stuff evrthg back
         n(i)=n1;
         r(i)=r1;
         r(i+1)=r2;
         t.n = n;
         t.r = r;
         t.core = cr;
         t.ps = psnew;
     otherwise
         error('Invalid subsasgn.');
 end


