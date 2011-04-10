function [elem] = subsref(tt,s)
%[ELEM]=SUBSREF(TT,S)
%Used to evaluate element of a tensor,
%and also --- get fields D,R,N,PS,CORE of the TT tensor

switch s(1).type    
    case '()'
	  %Evaluate element of a tensor in position s
      d=tt.d; n=tt.n; ps=tt.ps; cr=tt.core; r=tt.r;
      %ind=double(s);
      pp=s.subs;
      mn=numel(pp); 
      if ( mn < d && mn > 1 || (mn == 1 && numel(pp{1}) < d) ) 
        error('Invalid number of indices specified: given %d, need %d \n',mn,d);
      end
      if ( mn == 1 ) %Special case
         ind=pp{1};
         for i=1:d
           pp{i}=ind(i);
         end
      end
      elem=tt_tensor;
      elem.r=r; 
      elem.d=d;
      n0=zeros(d,1);
      for i=1:d
        ind=pp{i};
        if ( ind == ':')
            pp{i}=1:n(i);
            ind=pp{i};
        end
        n0(i)=numel(ind);
      end
      ps1=cumsum([1;n0.*r(1:d).*r(2:d+1)]);
      cr1=zeros(ps1(d+1)-1,1);
      for i=1:d
         ind=pp{i}; 
         core=cr(ps(i):ps(i+1)-1);
         core=reshape(core,[r(i),n(i),r(i+1)]);
         core=core(:,ind,:);
         cr1(ps1(i):ps1(i+1)-1)=core(:);          
      end
      if ( prod(n0) == 1 ) %single element
         v=1;
         for i=1:d
            core=cr1(ps1(i):ps1(i+1)-1); core=reshape(core,[r(i),r(i+1)]);
            v=v*core;
         end
         elem=v;
      else
         elem.n=n0;
         elem.core=cr1;
         elem.ps=ps1;
         
      end
    case '.'
          switch s(1).subs
             case 'r'
                 elem=tt.r;
             case 'core'
                 elem=tt.core;
             case 'd'
                 elem=tt.d;
              case 'ps' 
                  elem=tt.ps;
              case 'n'
                  elem=tt.n;
             otherwise
                 error(['No field ', s.subs, ' is here.']);
          end
    otherwise
        error('Invalid subsref.');
end
