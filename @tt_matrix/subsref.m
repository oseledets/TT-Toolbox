function elem = subsref(tt,s)
%[ELEM]=SUBSREF(TT,S)
%Used to evaluate element of a tensor

switch s(1).type
    case '()'
        error('Elementwise computation is not yet supported');
    case '.'
        switch s(1).subs               
            case 'n'
                elem=tt.n;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;                
            case 'm'
                elem=tt.m;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;    
            case 'tt'
                elem=tt.tt;
            otherwise
                % Dispatch to the underlying TT - may be rank, core, etc.
                elem = subsref(tt.tt, s);
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;                
%                 error(['No field ', s.subs, ' is here.']);
        end
    case '{}'
        pp=s.subs;
        mn=numel(pp);
        if ( mn > 1 )
          error('Invalid number of cores asked');
        end
        pp = pp{1};
        elem=core(tt.tt,pp);
        elem=reshape(elem, tt.tt.r(pp), tt.n(pp), tt.m(pp), tt.tt.r(pp+1));
        
        if (numel(s)>1)
            s = s(2:end);
            elem = subsref(elem, s);
        end;
    otherwise
        error('Invalid subsref.');
end
