function elem = subsref(tt,s)
% Evaluate cores of TT-matrix and fields of the TT-matrix structure
%   A=TT{I} computes the I-th core of the TT-representation
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
