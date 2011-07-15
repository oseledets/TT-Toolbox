function elem = subsref(tt,s)
%[ELEM]=SUBSREF(TT,S)
%Used to evaluate element of a tensor

switch s(1).type    
    case '()'
      error('Elementwise computation is not yet supported');
      case '.'
          switch s(1).subs
             case 'r'
                 elem=tt.tt.r;
             case 'tt'
                 elem=tt.tt;
             case 'd'
                 elem=tt.tt.d;
              case 'ps' 
                  elem=tt.tt.ps;
              case 'ns'
                  elem=tt.ns;
%               case 'm'
%                   elem=tt.m;
             otherwise
                 error(['No field ', s.subs, ' is here.']);
          end
    otherwise
        error('Invalid subsref.');
end
