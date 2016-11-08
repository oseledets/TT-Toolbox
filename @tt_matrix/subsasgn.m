function t = subsasgn(t,s,b)
%Indexing of cores of TT-format: T{I}=B
%   T{I} = B Set i-th core to be equal to B
%
%
%TT-Toolbox 2.2.2, 2009-2016
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%latest tweaks by Alexey Boyko,
%Skolkovo Institute of Science and Technology,
%a.boyko@skoltech.ru or alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2016
%---------------------------

switch s(1).type
    case '()'
        if size(s(1).subs)==[1 2] 
            if numel(s(1).subs{1})==1 && numel(s(1).subs{2})==1
                list_of_idx=1:numel(t.m);
                list_of_vals(1,:)=indexify(s(1).subs{1},t.n.');
                list_of_vals(2,:)=indexify(s(1).subs{2},t.m.');
            end
        end
        elem=tt_submatrix(t,list_of_idx,list_of_vals);
        t=t+(b-elem)*tt_onlyone([s(1).subs{1} s(1).subs{2}],t);
    case '.'
        switch s(1).subs
            case 'n'
                t.n=b;
            case 'm'
                t.m=b;
            case 'r'
                t.tt.r=b;
            case 'core'
                % I don't want to dispatch all the rest subs to the core,
                % because MATLAB will copy it. Fortunately, only () is
                % reasonable here.
                if (numel(s)>1)&&(strcmp(s(2).type,'()'))
                    if (size(s(2).subs{1},1)==1)
                        s(2).subs{1}=s(2).subs{1}';
                    end;
                    t.tt.core(s(2).subs{1})=b;
                else
                    t.tt.core=b;
                end;
            case 'ps'
                t.tt.ps=b;
            case 'tt'
                t.tt=b;
            case 'd'
                t.tt.d=b;
            otherwise
                error(['Cannot change field ', s.subs, ' directly.']);
        end
    case '{}'
        i = s.subs{1};
        r1 = size(b,1);
        n1 = size(b,2);
        m1 = size(b,3);
        r2 = size(b,4);
        b = reshape(b, r1, n1*m1, r2);
        t.tt{i} = b;
        t.n(i) = n1;
        t.m(i) = m1;
    otherwise
        error('Invalid subsasgn.');
end

end
