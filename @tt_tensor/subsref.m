function [elem] = subsref(tt,s)
% SUBSREF V2 
% Evaluate elements, cores of TT-formats and fields of the TT-structure
%   A=TT{I} computes I-th core of the TT-representation
%
%   ELEM=TT(IND) computes element with index IND of the tensor TT
%   
% TT-Toolbox 2.2.1, 2009-2016
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com

% Possibility to extract subtensor added by Alexey Boyko,
% Skolkovo Institute of Science and Technology
% For questions email alexey.boyko@skolkovotech.ru or a.boyko@skoltech.ru

% Example 1 
% y =tt_x(2,10)
% y2=y(1,2,1,2,2,:,:,:,:,:), or, alternatively, 
% y2=y(1,2,1,2,2) or
% y2=y([1,2,1,2,2]) or
% y2=y({1,2,1,2,2,':',':',':',':',':'}) or
% y2=y({1,2,1,2,2}) or
% will return 5-dimensional tt_tensor, corresponding to 
% fist 5 modes\harmonics being frozen to indices 1,2,1,2,2 respectively

% Example 2 
% y =tt_x(2,10)
% y2=y(1,:,1,2,2,:,:,1,:,:), or, alternatively, 
% y2=y({1,':',1,2,2,':',':',1}) or
% y2=y({1,':',1,2,2,':',':',1,':',':'}) or

% will also return 5-dimensional tt_tensor, corresponding to freezing
% respective modes\harmonics

%---------------------------

switch s(1).type    
    case '()'
    %Evaluate element of a tensor in position s
    d=tt.d; n=tt.n; ps=tt.ps; cr=tt.core; r=tt.r;

    %ind=double(s);
    pp=s.subs;
    mn=numel(pp); 

    if mn == 1
        if iscell( pp{1})
            %one argument, cell array
            pp=pp{1};
        else
            % one argument, ordinary array
            pp=num2cell(pp{1});
        end
    end
    mn=numel(pp); 

    
    
    if ( mn > d )
        error('Invalid number of indices specified: given %d, need no more than %d \n',mn,d);  
    end
    
    idxs=[];
    vals=[];
    for i=1:mn
        if isnumeric(pp{i})
            idxs=[idxs,i];
            vals=[vals,pp{i}];
        end
    end
    
    mn=numel(idxs);
    if(mn < d)
        elem=tt_subtensor(tt,idxs,vals);
    else

        %               if ( mn == 1 ) %Special case
        %                  ind=pp{1};
        %                  if ( is_array(ind) )
        %                    for i=1:d
        %                      pp{i}=ind(i);
        %                    end
        %                  else
        %                     pp=ind; %Supposedly it was a cell array of index subsets  
        %                  end
        %               end

        elem=tt_tensor;
        elem.r=r; 
        elem.d=d;
        n0=zeros(d,1);
        for i=1:d
            ind=pp{i};
            if ( strcmp(ind, ':'))
                pp{i}=1:n(i);
                ind=pp{i};
            end
            n0(i)=numel(ind);
        end
        ps1=cumsum([1;n0.*r(1:d).*r(2:d+1)]);
        cr1=zeros(ps1(d+1)-1,1);
        
        for i=1:d
            ind=pp{i}; 
            crx=cr(ps(i):ps(i+1)-1);
            crx=reshape(crx,[r(i),n(i),r(i+1)]);
            crx=crx(:,ind,:);
            cr1(ps1(i):ps1(i+1)-1)=crx(:);          
        end
        
        if ( prod(n0) == 1 ) %single element
            v=1;
            for i=1:d
                crx=cr1(ps1(i):ps1(i+1)-1); crx=reshape(crx,[r(i),r(i+1)]);
                v=v*crx;
            end
            elem=v;
        else
            elem.n=n0;
            elem.core=cr1;
            elem.ps=ps1;
        end
    
    end

    case '.'
        switch s(1).subs
            case 'r'
                elem = tt.r;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;
%                 if (numel(s)>1)&&(strcmp(s(2).type,'()'))
%                     elem = tt.r(s(2).subs{1});
%                 else
%                     elem=tt.r;
%                 end;
            case 'over'
                elem = tt.over;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;
            case 'core'
                if (numel(s)>1)&&(strcmp(s(2).type,'()'))
                    elem = tt.core(s(2).subs{1});
                else
                    elem=tt.core;
                end;
            case 'd'
                elem=tt.d;
            case 'ps'
                elem = tt.ps;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;

%                 if (numel(s)>1)&&(strcmp(s(2).type,'()'))
%                     elem = tt.ps(s(2).subs{1});
%                 else
%                     elem=tt.ps;
%                 end;
            case 'n'
                elem = tt.n;
                if (numel(s)>1)
                    s = s(2:end);
                    elem = subsref(elem, s);
                end;

%                 if (numel(s)>1)&&(strcmp(s(2).type,'()'))
%                     elem = tt.n(s(2).subs{1});
%                 else
%                     elem=tt.n;
%                 end;
            otherwise
                error(['No field ', s.subs, ' is here.']);
        end
        
    case '{}'
        %Return the core in the old (not exactly!!! r1-n-r2 here) format
        pp=s.subs;
        mn=numel(pp);
        if ( mn > 1 )
          error('Invalid number of cores asked');
        end
        elem=core(tt,pp{1});
        
        if (numel(s)>1)
            s = s(2:end);
            elem = subsref(elem, s);
        end;
%         if (pp{1}~=1)
%             elem=permute(elem,[2,1,3]);
%         end;
        
    otherwise
        error('Invalid subsref.');
end
end
