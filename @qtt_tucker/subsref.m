function [elem] = subsref(tt,s)
%Evaluate fields of the QTT-Tucker structure
%[ELEM]=SUBSREF(TT,S)
%Used to evaluate element of a tensor,
%and also --- get fields D,R,N,PS,CORE of the TT tensor

switch s(1).type    
    case '()'
        error('Element evaluation not implemented yet!');
    case '.'
        switch s(1).subs
            case 'core'
                elem = tt.core;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;                
            case 'tuck'
                elem = tt.tuck;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;
            case 'd'
                elem = tt.dphys;
            case 'dphys'
                elem = tt.dphys;
            case 'sz'
                elem = tt.sz;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;              
            otherwise
                error(['No field ', s.subs, ' is here.']);
        end
    case '{}'
%         %Return the core in the old (not exactly!!! r1-n-r2 here) format
        elem = subsref(tt.tuck, s);
%         pp=s.subs;
%         mn=numel(pp);
%         if ( mn > 1 )
%           error('Invalid number of cores asked');
%         end
%         elem=core(tt,pp{1});
% %         if (pp{1}~=1)
% %             elem=permute(elem,[2,1,3]);
% %         end;
        
    otherwise        
        error('Invalid subsref.');
end
