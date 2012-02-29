function [y]=reort(u,uadd)
%Golub-Kahan reorthogonalization
%   [Y]=REORT(U,UADD) Given orthonormal matrix U, orthonormalize UADD to it
%   to get orthogonal basis in [U,UADD]. Faster then the QR-decomposition,
%   using the Golub-Kahan method
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
if (size(uadd,2)==0)
    y = u;
    return;
end;
if (size(u,1) == size(u,2) )
  y=u;
  return
end
if (size(u,2) + size(uadd,2) >= size(u,1) )
  uadd=uadd(:,size(u,1)-size(u,2));
end

mvr=u'*uadd; unew=uadd-u*mvr; 
reort_flag=true;
while (reort_flag )
    %Nice vectorization!
    norm_unew=sum(unew.^2,1);
    norm_uadd=sum(uadd.^2,1);
    reort_flag=~isempty(find(norm_unew <= 0.25*norm_uadd,1));
    [unew,~]=qr(unew,0); %Here it is ok.
    if (reort_flag)
        su=u'*unew;
        uadd=unew;
        unew=unew-u*su;
    end
    
end
y=[u,unew];
return
end
