function [y]=funcrs2(tt,fun,eps,y,nswp)
%Cross approximation of a function of a TT-tensor, Method 2
%   [Y]=FUNCRS2(TT,FUN,EPS,Y,NSWP)
%   Computes approximation to the function FUN(TT) with accuracy EPS
%   Auxiliary parameters:  Y (initial approximation), NSWP
%   (number of sweeps in the method
%   Much faster then usual cross by vectorized computation of subtensors
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


%PARAMETERS SECTION
rmin=1;
verb=true;
kick_rank=5;
if (~isempty(y))
   yold=y;
end
y1=y;
%For this procedure we need only local (!) indices, since it
%is more or less equivalent to orthogonalization;
d=tt.d;
ps=tt.ps;
core0=tt.core;
n=tt.n;
r=tt.r;

ry=y.r;
psy=y.ps;
cry=y.core;
phi=cell(d+1,1);
phi{d+1}=1;
phi{1}=1;
phx=cell(d+1,1);
phx{d+1}=1; %For storing submatrices in U & V
phx{1}=1;
%Warmup procedure: orthogonalize from right-to-left & maxvol
   cry_old=cry; %This is for checking the accuracy
   psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
   pos1=psy(d+1);
   %The first is the the o
   for i=d-1:-1:1
      %Do right-to-left SVD + maxvol (i.e., no fun() is employed, just the
      %current approximation)
      cr=cry(pos1-ry(i+1)*n(i+1)*ry(i+2):pos1-1);
      cr2=cry_old(psy(i):psy(i+1)-1);
      cr2=reshape(cr2,[ry(i)*n(i),ry(i+1)]);
      cr=reshape(cr,[ry(i+1),n(i+1)*ry(i+2)]);
      %cr=cr.';
      %[cr,rm]=qr(cr,0);
      %rm=2*eye(ry(i+1));
      %cr=cr/2;
      [u,s,v]=svd(cr,'econ'); s=diag(s);
      ry(i+1) = my_chop2(s,norm(s)*eps/sqrt(d-1));
      %cr = u * s * v', I think we should leave v orthogonal
      u=u(:,1:ry(i+1)); v=v(:,1:ry(i+1));
      s=s(1:ry(i+1)); u=u*diag(s);
      cr=conj(v); %cr is n(i+1),ry(i+2),ry(i+1) ---  No. it is conj(v)
      rm=u.'; %This is discussable --- maybe u' (or u.')?
      ry(i+1)=size(cr,2);
      %Maxvol should be computed in a different matrix
      cr0=reshape(cr,[n(i+1),ry(i+2),ry(i+1)]);
      cr0=permute(cr0,[3,1,2]);
      cr0=reshape(cr0,[ry(i+1)*n(i+1),ry(i+2)]);
      cr0=cr0*phx{i+2}.';
      cr0=reshape(cr0,[ry(i+1),n(i+1)*ry(i+2)]);
      cr0=cr0.';
      ind=maxvol2(cr0);
      r1=cr0(ind,:);
      phx{i+1}=r1;

      cr=cr.';

      cry(pos1-ry(i+1)*n(i+1)*ry(i+2):pos1-1)=cr(:);
      pos1=pos1-ry(i+1)*n(i+1)*ry(i+2);
      cr2=cr2*(rm).';
      cry(pos1-ry(i)*n(i)*ry(i+1):pos1-1)=cr2(:);
      %Take phi matrix; convolve from right with current cors of V
      cr0=core0(ps(i+1):ps(i+2)-1);
      cr0=reshape(cr0,[r(i+1)*n(i+1),r(i+2)]);
      cr0=cr0*phi{i+2}; %cr0 is now r(i)*n(i)*ry(i+1);
      cr0=reshape(cr0,[r(i+1),n(i+1)*ry(i+2)]);
      phi{i+1}=cr0(:,ind);
   end
   pos1=pos1-ry(1)*n(1)*ry(2);
   psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
   cry=cry(pos1:numel(cry));
   y.core=cry;
   y.r=ry;
   y.ps=psy;

   swp=1;
yold=[];
not_converged = true;
while ( swp < nswp && not_converged )
    max_er=0;
cry_old=cry; %This is for checking the accuracy
psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
pos1=1;
 for i=1:d-1
     %fprintf('i=%d \n',i);
     %We care for two cores, with number i & number i+1, and use
     %psi(i) and psi(i+2) as a basis; also we will need to recompute
     %psi(i+1)
     ps1=phi{i}; ps2=phi{i+2}; px1=phx{i}; px2=phx{i+2};
     %Compute (!) superblock and function (!) of it
     cr1=core0(ps(i):ps(i+1)-1);
     cr2=core0(ps(i+1):ps(i+2)-1);
     cr1=reshape(cr1,[r(i),n(i)*r(i+1)]);
     cr1=ps1*cr1;
     cr2=reshape(cr2,[r(i+1)*n(i+1),r(i+2)]);
     cr2=cr2*ps2;
     cr1=reshape(cr1,[ry(i)*n(i),r(i+1)]);
     cr2=reshape(cr2,[r(i+1),n(i+1)*ry(i+2)]);
     cr=cr1*cr2;
     %cr=reshape(cr,[ry(i)*n(i)*n(i+1),ry(i+2)]);
     cr=fun(cr); %Elements are evaluated here!
     cr=reshape(cr,[ry(i)*n(i)*n(i+1),ry(i+2)]);
     cr = cr/(px2.');
     cr=reshape(cr,[ry(i),n(i)*n(i+1)*ry(i+2)]);
     cr=px1 \cr;
     cr=reshape(cr,[ry(i)*n(i),n(i+1)*ry(i+2)]);
     %Check for local approximation of cr for the error
     cry1=cry(pos1:pos1+ry(i)*n(i)*ry(i+1)-1);
     %cry2=cry(pos1+ry(i)*n(i)*ry(i+1):pos1+ry(i)*n(i)*ry(i+1)+ry(i+1)*n(i+1)*ry(i+2)-1);
     %if ( swp == 2 )
     % if ( i == 5 )
     %   keyboard;
     % end
     %end
     cry2=cry_old(psy(i+1):psy(i+2)-1);

     cry1=reshape(cry1,[ry(i)*n(i),ry(i+1)]);
     cry2=reshape(cry2,[ry(i+1),n(i+1)*ry(i+2)]);
     appr=cry1*cry2;
     er=norm(appr-cr,'fro')/norm(cr,'fro');
     %keyboard;
     max_er=max(er,max_er);

     %Compute SVD of cr

     [u,s,v]=svd(cr,'econ');
     s=diag(s);
     r2=my_chop2(s,eps*norm(s)/sqrt(d-1));
     s=s(1:r2); u=u(:,1:r2); v=conj(v(:,1:r2));
     v=v*diag(s);

     %Kick rank of u
     ur=randn(size(u,1),kick_rank);
     %Orthogonalize ur to u by Golub-Kahan reorth <- it sucks!!
     [u,rv]=qr([u,ur], 0);
%      u=reort(u,ur);
%      radd=size(u,2)-r2;
%      if ( radd > 0 )
       vr=zeros(size(v,1),kick_rank);
       v=[v,vr];
       v = v*(rv.');
%      end
     r2=size(u,2);


     %ind=maxvol2(u);
     %r1=u(ind,:);
     %u=u/r1; v=v*r1'; v=v.';
     v=v.';
     ry(i+1)=r2;
     cry(pos1:pos1+ry(i)*n(i)*ry(i+1)-1)=u(:);
     pos1=pos1+ry(i)*n(i)*ry(i+1);
     cry(pos1:pos1+ry(i+1)*n(i+1)*ry(i+2)-1)=v(:);



     %Compute new maxvol and (new) submatrices
     u=reshape(u,[ry(i),n(i)*ry(i+1)]);
     u=px1*u;
     u=reshape(u,[ry(i)*n(i),ry(i+1)]);
     ind=maxvol2(u);
     phx{i+1}=u(ind,:);
     %Now we have to: cr1 with phi from the left
     phi{i+1}=cr1(ind,:); %phi{i+1}=phi{i+1};


 end
  %Truncate local memory
  cry=cry(1:pos1+ry(d)*n(d)*ry(d+1)-1);
  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
y.core=cry;
y.r=ry;
y.ps=psy;
cry_old=cry;
%m=numel(cry);
%cry=[zeros(size(cry)),cry];
%pos1=pos1+m;
%cry2=cry_old(psy(i+1):psy(i+2)-1);
cry=cry(psy(d):psy(d+1)-1); %Start--only two cores
 for i=d-1:-1:1
     %fprintf('i=%d \n',i);
     %We care for two cores, with number i & number i+1, and use
     %psi(i) and psi(i+2) as a basis; also we will need to recompute
     %psi(i+1)
     ps1=phi{i}; ps2=phi{i+2}; px1=phx{i}; px2=phx{i+2};
     %Take current core; convolve it with
     %core=core(ps
     %Compute (!) superblock and function (!) of it
     cr1=core0(ps(i):ps(i+1)-1);
     cr2=core0(ps(i+1):ps(i+2)-1);
     cr1=reshape(cr1,[r(i),n(i)*r(i+1)]);
     cr1=ps1*cr1;
     cr2=reshape(cr2,[r(i+1)*n(i+1),r(i+2)]);
     cr2=cr2*ps2;
     cr1=reshape(cr1,[ry(i)*n(i),r(i+1)]);
     cr2=reshape(cr2,[r(i+1),n(i+1)*ry(i+2)]);
     cr=cr1*cr2;
     cr=reshape(cr,[ry(i)*n(i),n(i+1)*ry(i+2)]);
     cr=fun(cr); %Elements are evaluated here!
     cr=reshape(cr,[ry(i)*n(i)*n(i+1),ry(i+2)]);
     cr = cr/(px2.');
     cr=reshape(cr,[ry(i),n(i)*n(i+1)*ry(i+2)]);
     cr=px1 \cr;
     cr=reshape(cr,[ry(i)*n(i),n(i+1)*ry(i+2)]);
     %Check for local approximation of cr for the error
     %cry1=cry(pos1-ry(i)*n(i)*ry(i+1):pos1-1);
     cry1=cry_old(psy(i):psy(i+1)-1);
     %cry2=cry(ry(:ry(i+1)*n(i+1)*ry(i+2));
     %cry1=cry(1:ry(i)*n(i)*ry(i+1));
     %cry2=cry(ry(i)*n(i)*ry(i+1)+1:ry(i)*n(i)*ry(i+1)+ry(i+1)*n(i+1)*ry(i+2));
     cry2=cry(1:ry(i+1)*n(i+1)*ry(i+2));
     cry(1:ry(i+1)*n(i+1)*ry(i+2))=[];
     %cry2=cry_old(psy(i+1):psy(i+2)-1);
     %cry(1:(ry(i+1)*n(i+1)*ry(i+2)))=[]; %Delete both of first cores
     %cry(1:ry(i)*n(i)*ry(i+1)+ry(i+1)*n(i+1)*ry(i+2))=[];
     cry1=reshape(cry1,[ry(i)*n(i),ry(i+1)]);
     cry2=reshape(cry2,[ry(i+1),n(i+1)*ry(i+2)]);
     appr=cry1*cry2;
     er=norm(appr-cr,'fro')/norm(cr,'fro');
     %er
          max_er=max(er,max_er);

     %Compute SVD of cr
     [u,s,v]=svd(cr,'econ');
     s=diag(s);
     r2=my_chop2(s,eps*norm(s)/sqrt(d-1));
     s=s(1:r2); u=u(:,1:r2); v=v(:,1:r2);
     
     v = conj(v); % <- bug was here
     
     %Make it standard
     u=u*diag(s);


     %Kick rank

      vr=randn(size(v,1),kick_rank);
      [v,rv]=qr([v,vr], 0);
%       v=reort(v,vr);
%       radd=size(v,2)-r2;
%       if ( radd > 0 )
         ur=zeros(size(u,1),kick_rank);
         u=[u,ur];
         u = u*(rv.');
%       end
      r2=size(v,2);



     v=v.';v0=v;
     %Compute new phi;
     ry(i+1)=r2;
     u=u(:); u=u.'; v=v(:); v=v.';
     cry=[u,v,cry];
     %keyboard;
     %We need new memory for
     %cry(pos1:pos1+ry(i+1)*n(i+1)*ry(i+2)-1)=v(:); %Here memory has to be (?) enlarged
     %if ( pos1 <= ry(i)*n(i)*ry(i+1) ) %Have to enlarge memory
     %  cry=cry(pos1:numel(cry));
     %  cry=[u(:)',cry];
     %  pos1=ry(i)*n(i)*ry(i+1)+1;
     %else
     %   pos1=pos1-ry(i)*n(i)*ry(i+1);
      %  cry(pos1:pos1+ry(i)*n(i)*ry(i+1)-1)=u(:);
     %end
     %Compute new maxvol and (new) submatrices
     v0=reshape(v0,[ry(i+1)*n(i+1),ry(i+2)]);
     v0=v0*px2.';
     v0=reshape(v0,[ry(i+1),n(i+1)*ry(i+2)]);
     v0=v0.';
     ind=maxvol2(v0);
     phx{i+1}=v0(ind,:);
     %Now we have to: cr1 with phi from the left
     phi{i+1}=cr2(:,ind);

 end

  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);

y.core=cry;
y.r=ry;
y.ps=psy;
if ( isempty(yold) )
 yold=y;
 er_nrm=1;
else

   er_nrm=norm(yold-y)/norm(y);
   yold=y;
end
 fprintf('sweep=%d, er=%3.2e er_nrm=%3.2e \n',swp,max_er,er_nrm);

swp=swp+1;
end
  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);

y.core=cry;
y.r=ry;
y.ps=psy;

end
