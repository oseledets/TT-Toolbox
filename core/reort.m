function [y]=reort(u,uadd)
%[Y]=REORT(U,UADD)
%Faster (?) QR-decomposition of the matrix [u,v] by Golub-Kahan
%reorthogonalization
if (size(u,1) == size(u,2) )
  y=u;
  return
end
if (size(u,2) + size(uadd,2) >= size(u,1) )
  %y=eye(size(u,1)); %The simplest basis
  uadd=uadd(:,size(u,1)-size(u,2));
end

mvr=u'*uadd; unew=uadd-u*mvr; 
reort_flag=true;
radd=size(unew,2);
while (reort_flag )
    reort_flag=false;
for j=1:radd
  if ( norm(unew(:,j)) <= 0.5*norm(uadd(:,j)))
          reort_flag=true;
  end
end
if (reort_flag)
  su=u'*unew;
  uadd=unew;
  unew=unew-u*su; 
end
  
end
[unew,~]=qr(unew,0); 
y=[u,unew];

return
end