function t = subsasgn(t,s,b)
%Assign fields of the QTT-Tucker structure
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
     otherwise
         error('Invalid subsasgn.');
 end


