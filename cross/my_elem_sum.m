function [val]=my_elem(ind,arr)
%Prepare for TT-cross method applied to symplex
%true_d=arr{1}
%true_size=arr{2};
K   = arr{3};
d   = arr{4}; 
h   = arr{5};
a   = arr{6};
prm = arr{7};

%ind=ind(prm);
x=zeros(K,1);
for k=1:K
    indk=ind((k-1)*d+1:d*k);
  indk = sub_to_ind(indk,2*ones(1,d));
  x(k) = a + indk*h;

end
%val=1.0/(norm(x)+1e-6);
val=sum(x);
%val=sum(x);
%if ( val <= 1 )
%  val=1.0;
%else
%  val=0.0;
%end
return
end
