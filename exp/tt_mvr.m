%%
function [y,phi] = tt_mvr(a,x,y)
% [RES] = TT_MVR(A,X,NSWP)
% This is a technical function, tbat computes initialization for
% TT_MVK
%
% TT Toolbox 1.0, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
d=size(a,1);
phi=cell(d,1); 
%phi is rx x ra
%phi will be matrices
%The first step is to compute left-to-right sweep of summation, yielding
%phi matrices
for i=1:d
  y{i}=y{i}./norm(y{i});
end
%A(n,m,ra)*X(m,rx)*Y(n)
%Decipher sizes
core=a{1}; n=size(core,1); m=size(core,2); ra=size(core,3);
matx=x{1}; rx=size(matx,2);
vecy=y{1};
ph=reshape(core,[n,m*ra])'*vecy;
%ph is m*ra
ph=matx'*reshape(ph,[m,ra]);
phi{1}=ph;
for i=2:d-1
   %phi{i-1} is rx1 x ra1
   %phi1(rx1,ra1)*A(n,m,ra1,ra2)*X(m,rx1,rx2)*Y(n)
   core=a{i}; core_x=x{i}; vecy=y{i};
   %Decipher size
   n=size(core,1); m=size(core,2); ra1=size(core,3); ra2=size(core,4);
   rx1=size(core_x,2); rx2=size(core_x,3);
   ph=vecy'*reshape(core,[n,m*ra1*ra2]); %ph is m*ra1*ra2
   %phi(m,ra1,ra2)*phi1(rx1,ra1)
   ph=reshape(ph,[m,ra1,ra2]); ph=permute(ph,[2,1,3]);
   ph=phi{i-1}*reshape(ph,[ra1,m*ra2]); %ph is rx1 x m x ra2
   ph=reshape(ph,[rx1*m,ra2]);
  
   %end
   ph=reshape(permute(core_x,[2,1,3]),[rx1*m,rx2])'*ph; %ph rx2*ra2
   phi{i}=ph;
end
%Now we are able to recompute the last vector
%A(n,m,ra)*X(m,rx)*phi(rx,ra) -> new y (and new sigma!)
   %Decipher sizes
   core=a{d}; matx=x{d};
   ph=phi{d-1}; 
   n=size(core,1); m=size(core,2); ra=size(core,3); rx=size(matx,2);
   ph=matx*ph;
   vecy=reshape(core,[n,m*ra])*reshape(ph,[m*ra,1]); 
   sigma=norm(vecy);
   y{d}=vecy./sigma;
   %And recompute phi:
   %A(n,m,ra)*X(m,rx)*Y(n)=phi(rx,ra)
   core=a{d}; n=size(core,1); m=size(core,2); ra=size(core,3);
   matx=x{d}; rx=size(matx,2);
   yvec=y{d};
   ph=reshape(core,[n,m*ra])'*yvec;
%ph is m*ra
   ph=matx'*reshape(ph,[m,ra]);
   phi{d}=ph;
   for i=d-1:-1:2
       %We compute
       %phi1(rx1,ra1)*A(n,m,ra1,ra2)*X(m,rx1,rx2)*phi2(rx2,ra2)
       %Plus A(n,m,ra1,ra2)*X(m,rx1,rx2)*phi(rx2,ra2)*Y(n) -> (rx1,ra1)
       %Therefore, first compute
       %A(n,m,ra1,ra2)*X(m,rx1,rx2)*phi2(rx2,ra2)
       core=a{i}; core_x=x{i}; ph2=phi{i+1}; ph1=phi{i-1};
       n=size(core,1); m=size(core,2); ra1=size(core,3); ra2=size(core,4);
       rx1=size(core_x,2); rx2=size(core_x,3); 
       ph2=reshape(core_x,[m*rx1,rx2])*ph2; %ph2 is m*rx1*ra2
       ph2=reshape(ph2,[m,rx1,ra2]); ph2=permute(ph2,[2,1,3]); 
       core=permute(core,[1,3,2,4]); core=reshape(core,[n*ra1,m*ra2]);
       ph2=core*reshape(ph2,[rx1,m*ra2])';
       %ph2=reshape(permute(core,[1,3,2,4]),[n*ra1,m*ra2])*reshape(ph2,[rx1,m*ra2])';
       %ph2 is n*ra1*rx1 -> this is required
       %Compute  new vector
       ph2=reshape(ph2,[n,ra1,rx1]); ph2=permute(ph2,[1,3,2]); %here it is n*rx1*ra1
       vecy=reshape(ph2,[n,rx1*ra1])*reshape(ph1,[rx1*ra1,1]);
       sigma=norm(vecy);
       y{i}=vecy./sigma;
       %And new phi
       ph1=y{i}'*reshape(ph2,[n,rx1*ra1]);
       phi{i}=reshape(ph1,[rx1,ra1]);       
   end
   %Finally, compute y{1};
   %A(n,m,ra)*X(m,rx)*phi(rx,ra)
   %and 
   %A(n,m,ra)*X(m,rx)*Y(n)
   core=a{1}; matx=x{1}; ph=phi{2};
   n=size(core,1); m=size(core,2); ra=size(core,3);
   rx=size(matx,2);
   ph=matx*ph; %is m*ra
   vecy=reshape(core,[n,m*ra])*reshape(ph,[m*ra,1]);
   sigma=norm(vecy);
   y{1}=vecy./sigma;
   %vecy=y{1};
   %
   %ph=reshape(core,[n,m*ra])'*vecy;
   %ph is m*ra
   %ph=matx'*reshape(ph,[m,ra]);
   %phi{1}=ph;
   
return;
end
