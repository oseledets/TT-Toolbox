function t = subsasgn(t,s,b)
% 
 switch s(1).type
     case '.'
         switch s(1).subs
             case 'dphys'
                 t.dphys = b;
             case 'd'
                 t.dphys = b;                 
             case 'core'
                 if (numel(s)>1)
                     % Dispatch extra indexing gracefully
                     t.core = subsasgn(t.core, s(2:end), b);
                 else
                     t.core = b;
                 end;                 
             case 'tuck'
                 if (numel(s)>1)
                     % Dispatch extra indexing gracefully <-BULLSHIT! IT
                     % DOESN'T WORK!!!
                     s = s(2:end);
                     if (strcmp(s(1).type, '{}'))
                         i = s(1).subs{1};
                         if (numel(s)>1)
                             t.tuck{i} = subsasgn(t.tuck{i}, s(2:end), b);
                         else
                             t.tuck{i} = b;
                         end;
                     end;
                 else
                     t.tuck = b;
                 end;
             case 'sz'
                 if (numel(s)>1)
                     % Dispatch extra indexing gracefully
                     t.sz = subsasgn(t.sz, s(2:end), b);
                 else
                     t.sz = b;
                 end;
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


