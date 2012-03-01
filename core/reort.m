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
  uadd=uadd(:,1:(size(u,1)-size(u,2)));
end

mvr=u'*uadd; unew=uadd-u*mvr; 
reort_flag=true;
j=1;
while (reort_flag && j <= 2)
    %Nice vectorization!
    norm_unew=sum(unew.^2,1);
    norm_uadd=sum(uadd.^2,1);
    %if (sum(norm_unew,2)<1e-28*sum(norm_uadd,2))
    %    unew = [];
    %    break;
    %end;
    reort_flag=~isempty(find(norm_unew <= 0.25*norm_uadd,1));
    [unew,~]=qr(unew,0); %Here it is ok.
    if (reort_flag)
        su=u'*unew;
        unew=unew-u*su;
        %Kill machine zeros
        unew(unew==0)=max(1e-308,norm(unew,'fro')*1e-15);
        j=j+1;
    end
end
if ( reort_flag )
 fprintf('Reort failed to reort! \n');
 [y,~]=qr([u,unew]);
else
 y=[u,unew];
%end
return
end
