function elem = subsref(tt,s)
% Evaluate cores of TT-matrix and fields of the TT-matrix structure
%   A=TT{I} computes the I-th core of the TT-representation
%   
%
% TT-Toolbox 2.2.2, 2009-2016
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
        if size(s(1).subs)==[1 2] 
            if numel(s(1).subs{1})==1 && numel(s(1).subs{2})==1
                list_of_idx=1:numel(tt.m);
                list_of_vals(1,:)=indexify(s(1).subs{1},tt.n.');
                list_of_vals(2,:)=indexify(s(1).subs{2},tt.m.');
            elseif size( s(1).subs{2},1 )==2 && size( s(1).subs{1},2 )==size( s(1).subs{2},2 )
                if size( s(1).subs{1},1 )==1
                    list_of_idx=s(1).subs{1};
                    list_of_vals=s(1).subs{2};
                elseif size( s(1).subs{1},1 )==2
                    error('Assymetric submatrix extraction is not yet supported');
                end
            end
        end
        elem=tt_submatrix(tt,list_of_idx,list_of_vals);
%         %ind=double(s);
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
