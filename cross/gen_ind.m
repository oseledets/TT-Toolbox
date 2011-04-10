%Generates index sets for the TT decomposition
function [ind]=gen_ind(sz,opt,r)
if ( size(sz,1) == 1 )
  d = size(sz,2);
else
  d = size(sz,1);
end
ind=cell(d,1);
if ( opt == 'r' )
    %From the right
    for k=1:d-1
        q=prod(sz((k+1):d));
        q=min(q,r);
        ind_cur=zeros(d-k,q);
        for i=1:q
            %ind_f=i;
            ind_cur(:,i)=ones(d-k,1);%ind2sub(sz((k+1):d),ind_f);
        end
        ind{k}=ind_cur;
    end
elseif ( opt == 'l')
else
    fprintf('gen_ind: Incorrect option \n');
    return
end
    
return
end