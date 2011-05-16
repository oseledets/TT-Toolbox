function [y]=funcross(tt,fun,eps,y,nswp)
%[Y]=FUNCROSS(TT,FUN,EPS)
%[Y]=FUNCROSS(TT,FUN,EPS)
%[Y]=FUNCROSS(TT,FUN,EPS,Y)
%[Y]=FUNCROSS(TT,FUN,EPS,Y,NSWP)
%Computes approximation to the function FUN(TT) with accuracy EPS
%Auxiliary parameters:  Y (initial approximation), NSWP 
%(number of sweeps in the method 
%Much faster then usual cross by vectorized computation of subtensors
%It is now based on the "add" concept, not on the "replace" concept 
%for the supercore. The compression is performed in a backwards sweep.

%The scheme of the algorithm:
%(0) --- warmup step, l-r qr. of the initial guess (no psi are needed)
%MAIN CYCLE
%Compress from right-to-left & compute maxvol + psi (?) matrices
%Compute left-to-right DMRG maxvol in the add scheme,
%The idea is simple: each supercore U*V' is randomly pushed as
%[U,u]*[V,0]'





%PARAMETERS SECTION
verb=true;
ry_old=5; %This is the number of random vectors to be added    
if (~isempty(y))
   yold=y;
end



%For this procedure we need only local (!) indices, since it
%is more or less equivalent to orthogonalization; 
d=tt.d; 
ps=tt.ps;
core0=tt.core;
n=tt.n;
r=tt.r;
swp=1;
phi=cell(d+1,1);
phi{d+1}=1;
phi{1}=1;
phx=cell(d+1,1);
phx{d+1}=1; %For storing submatrices in U & V
phx{1}=1;
%First step is to orthogonalize the input from left-to-right;
%this is the current "warm up" step;
y=qr(y,'lr');
ry=y.r;
cry=y.core;

not_converged = true;
converged_once=false;
while ( swp < nswp && not_converged )
   %y0=tt_random(size(y),y.d,y.r(2:y.d));
   %y0=tt_random(size(y),y.d,2);
   
   %y0=tt_tensor(y0); 
   %y=y+y0;
   cry=y.core;
   ry=y.r;
    psy=y.ps;
   cry_old=cry; %This is for checking the accuracy
   pos1=psy(d+1);
   %The first is the the orthogonalization
   for i=d-1:-1:1
      %Do right-to-left SVD + maxvol (i.e., no fun() is employed, just the
      %current approximation)
      cr=cry(pos1-ry(i+1)*n(i+1)*ry(i+2):pos1-1);
      cr2=cry_old(psy(i):psy(i+1)-1);
      cr2=reshape(cr2,[ry(i)*n(i),ry(i+1)]);
      cr=reshape(cr,[ry(i+1),n(i+1)*ry(i+2)]);
      
      
      [u,s,v]=svd(cr,'econ'); s=diag(s); 
       %ry(i+1) = my_chop2(s,norm(s)*eps/sqrt(d-1));
       ry(i+1)=my_chop2(s,-1);
       ry(i+1)=min(ry(i+1),numel(s));
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
   %Now left-to-right (and double the core)
  pos1=1;
  max_er=0;
  cry_old=cry; %This is for checking the accuracy
 for i=1:d-1
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
     cry2=cry_old(psy(i+1):psy(i+2)-1);
     
     cry1=reshape(cry1,[ry(i)*n(i),ry(i+1)]);
     cry2=reshape(cry2,[ry(i+1),n(i+1)*ry(i+2)]);
     
     appr=cry1*cry2;
     er=norm(appr-cr,'fro')/norm(cr,'fro');

     max_er=max(er,max_er);
     
     %This is the random add step. Serves to improve the convergence
     
     %if ( converged_once )
     %  ry_add=ry(i+1);
     %else
       ry_add=ry_old;
     %end
     cry1=randn(ry(i),n(i),ry_add); 
     cry1=reshape(cry1,[ry(i)*n(i),ry_add]);
     [cry1,~]=qr(cry1,0);
     ry_add=size(cry1,2);
     %ttt=randn([ry(i),n(i),ry(i+1)+ry_old]);
     %ttt(:,:,1:ry(i+1))=reshape(cry1,ry(i),n(i),ry(i+1));
     %cry1=ttt;
     %ry_add=ry(i+1)+ry_old;
     
     
     
       
     %Compute SVD of cr

     [u,s,v]=svd(cr,'econ'); 
     s=diag(s);
     r2=my_chop2(s,eps*norm(s)/sqrt(d-1));
     s=s(1:r2); u=u(:,1:r2); v=v(:,1:r2);
     v=v*diag(s);

     v=v.';
     ry(i+1)=r2;
    
     
     unew=zeros(ry(i),n(i),ry(i+1)+ry_add);
     vnew=zeros(ry_add+ry(i+1),n(i+1),ry(i+2));
     unew(:,:,1:ry(i+1))=reshape(u,[ry(i),n(i),ry(i+1)]);
     unew(:,:,ry(i+1)+1:end)=reshape(cry1,[ry(i),n(i),ry_add]);
      
     vnew(1:ry(i+1),:,:)=reshape(v,[ry(i+1),n(i+1),ry(i+2)]);
     
     %unew=u;
     %vnew=v;
     
     %Orthogonalize unew
     ry(i+1)=ry(i+1)+ry_add;
     [unew,rm]=qr(reshape(unew,[ry(i)*n(i),ry(i+1)]),0);
     vnew=rm*reshape(vnew,[ry(i+1),n(i+1)*ry(i+2)]);
     ry(i+1)=size(unew,2);
     cry(pos1:pos1+ry(i)*n(i)*ry(i+1)-1)=unew(:);
     pos1=pos1+ry(i)*n(i)*ry(i+1);
     cry(pos1:pos1+ry(i+1)*n(i+1)*ry(i+2)-1)=vnew(:);
     
    
     
     %Compute new maxvol and (new) submatrices
     unew=reshape(unew,[ry(i),n(i)*ry(i+1)]);
     unew=px1*unew;
     unew=reshape(unew,[ry(i)*n(i),ry(i+1)]);
     ind=maxvol2(unew);
     phx{i+1}=unew(ind,:);
     %Now we have to: cr1 with phi from the left 
     phi{i+1}=cr1(ind,:); %phi{i+1}=phi{i+1};
     
     
 end
  %Truncate local memory
  cry=cry(1:pos1+ry(d)*n(d)*ry(d+1)-1);
  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
y.core=cry;
y.r=ry;
y.ps=psy;
if ( isempty(yold) )
 yold=y;
 er_nrm=1;
else
   
   er_nrm=norm(yold-y)/norm(y);
   if ( isinf(er_nrm) || isnan(er_nrm) )
  %    keyboard;
   end
   yold=y;
end
if ( verb )
    
 fprintf('sweep=%d, er=%3.2e er_nrm=%3.2e \n',swp,max_er,er_nrm);
end
if ( max_er < eps && er_nrm < eps )
    if ( converged_once )
      not_converged=false;
    else
      converged_once=true;    
      ry_old=ry_old*2;
    end
else
  converged_once=false;
end
swp=swp+1;
end     
  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);

y.core=cry;
y.r=ry;
y.ps=psy;

end
