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
                t.n=b;
            case 'm'
                t.m=b;
            case 'r'
                t.tt.r=b;
            case 'core'
                t.tt.core=b;
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