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
%Scale uadd
radd=size(uadd,2);

mvr=u'*uadd; unew=uadd-u*mvr; 
reort_flag=true;
while (reort_flag )
    reort_flag=false;
    j=1;
while ( ~reort_flag && j <= radd)
  if ( norm(unew(:,j)) <= 0.5*norm(uadd(:,j)))
          reort_flag=true;
  end
  j=j+1;  
end
[unew,~]=qr(unew,0); %Here it is ok.
if (reort_flag)
  su=u'*unew;
  uadd=unew;
  unew=unew-u*su; 
  %[unew,~]=qr(unew,0);
end
  
end
y=[u,unew];
%if ( norm(y'*y-eye(size(y,2))) > 1e-13 )
%  keyboard
%end
return
end