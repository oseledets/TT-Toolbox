function t = subsasgn(t,s,b)
% 
 switch s(1).type
     case '.'
         switch s(1).subs
             case 'n'
                t.n=b;
             case 'r'
                t.r=b;
             case 'core'
                t.core=b;
             case 'd'
                 t.d=b;
             case 'ps'
                 t.ps=b;
             otherwise
                 error(['Cannot change field ', s.subs, ' directly.']);
         end
     otherwise
         error('Invalid subsasgn.');
 end


