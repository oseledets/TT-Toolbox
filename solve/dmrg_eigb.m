function [y,ev] = dmrg_eigb(a,k,eps,y0,rmax,nswp)
%[y,ev]=DMRG_EIGB(A,K,EPS,[Y0],[RMAX],[NSWP])
%Solves for K minimal eigenvalues of the TT-matrix A 
%with accuracy EPS by minimizing block Rayleigh quotient
%by the two-sided DMRG method
%Y0 is a possible initial guess, RMAX is the maximal rank
%NSWP is the maximal number of sweeps
%Currently generalized eigenproblem is not supported (hope to 
%do it in future)
n=a.n; %This works only for square matrices :)
d=a.d;
if ( nargin <= 3 || isempty(y0) )
    %Generate a random block tensor with the block dimension on the
    %right 
   y=tt_tensor;
   y.d=d;
   ry=2*ones(1,d); ry(1)=1; ry(d+1)=k; ry=ry';
   y.r=ry;
   y.n=n;
   y.ps=cumsum([1;y.n.*ry(1:d).*ry(2:d+1)]);
   psy=y.ps;
   sz=psy(d+1)-1;
   cr=randn(1,sz);
   %cr=ones(sz,1);
   y.core=cr;
else
   y=y0;
end
if ( nargin <= 4 || isempty(rmax) )
  rmax=1000;
end
if ( nargin <= 5 || isempty(nswp) )
  nswp = 40;
end
y0=y;
%            fm=full(a); 
%            [v,dg]=eig(fm);
%            ev=diag(dg);
%            [ev,ind]=sort(ev,'ascend');
%            v=v(:,ind);
%           w=v(:,1:k); ww0=w;
%            ev=ev(1:k);
%  w=reshape(w,[2*ones(1,d),k]);  
%  y=tt_tensor(w,1e-8);
%  y.d=d;
%  psy=y.ps; 
%  psy=psy(1:d+1);
%  cry=y.core; 
%  cry=cry(1:psy(d+1)-1); 
%  y.core=cry;
%  y.n=2*ones(d,1);
%  ry=y.r; ry=ry(1:d+1); 
%  y.r=ry;


%Parameters section
msize=10;
max_l_steps=1000;
kick_rank=5;

%We start from the orthogonalization of the y vector from left-to-right
%(it does not influence the TT-ranks)

psy=y.ps;
ry=y.r;
tta=a.tt;
psa=tta.ps;
ra=tta.r;
cra=tta.core;
cry=y.core;
rm=1; 
%Also we will need to compute phi-matrices
phi=cell(d+1,1);
phi{1}=1; %This one is right now unchanged
phi{d+1}=1; %Seems to be unchanged also
pos1=1;
for i=1:d
    cr=cry(psy(i):psy(i+1)-1);
    cr=reshape(cr,[ry(i),n(i)*ry(i+1)]);
    cr=rm*cr; 
    ry(i)=size(cr,1); 
    cr=reshape(cr,[ry(i)*n(i),ry(i+1)]);
    [u,rm]=qr(cr,0);
    rnew=size(u,2); 
    %ry(i+1)=rnew;
    cry(pos1:pos1+ry(i)*n(i)*rnew-1)=u(:); 
    pos1=pos1+ry(i)*n(i)*rnew;
    %With this new core compute (Ax,x) to the phi matrix
    %u=reshape(u,[ry(i),n(i),ry(i+1)]);
    crm=cra(psa(i):psa(i+1)-1);
    crm=reshape(crm,[ra(i),n(i),n(i),ra(i+1)]); 
    phx=phi{i};
    phx=reshape(phx,[ra(i),ry(i),ry(i)]); phx=permute(phx,[1,3,2]);
    x0=u;
    x0=reshape(x0,[ry(i),n(i)*rnew]);
    phx=reshape(phx,[ra(i)*ry(i),ry(i)]);
    phx=phx*x0; %phx is ra(i)*ry(i)*n(i)*ry(i+1)
    phx=reshape(phx,[ra(i),ry(i),n(i),rnew]);
    crm=permute(crm,[1,3,2,4]); crm=reshape(crm,[ra(i)*n(i),n(i)*ra(i+1)]);
    %Convolve over (ak-1) j with the matrix
    phx=permute(phx,[2,4,1,3]); 
    phx=reshape(phx,[ry(i)*rnew,ra(i)*n(i)]);
    phx=phx*crm; %ry(i)*ry(i+1)*n(i)*ra(i+1)
    phx=reshape(phx,[ry(i),rnew,n(i),ra(i+1)]);
    phx=permute(phx,[2,4,1,3]);
    phx=reshape(phx,[rnew*ra(i+1),ry(i)*n(i)]);
    x0=u;
    x0=reshape(x0,[ry(i)*n(i),rnew]);
    phx=phx*x0; 
    phx=reshape(phx,[rnew,ra(i+1),rnew]); 
    phx=permute(phx,[2,1,3]);
    phi{i+1}=phx;
    %phi{i+1}=permute(phi{i+1},[1,3,2]);
end

%phh0=phi;
ry(d+1)=rnew;
 %Test back transform
%  y0=y;
%  y0.n=[2*ones(d,1);k];
%  y0.r=[y.r;1];
%  y0.d=d;
%  e0=eye(ry(d));
%  cry=[cry,e0(:)'];
%  y0.core=cry;
%  ww1=full(y0); ww1=reshape(ww1,[numel(ww1)/k,k]);
%  bm=ww0'*ww1;
%  norm(ww1-ww0*bm)
%keyboard;
phi{d+1}=1; %Bydlocode, but seems necessary
%Truncate the core
%pos1=pos1-ry(d)*n(d)*ry(d+1);
cry=cry(1:pos1-1);
y.r=ry;
y.core=cry;
y.ps=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
psy=y.ps;
%Now start the main DMRG sweep
swp=1;
not_converged=true;
%ry(1)=1;
ry(d+1)=1;
dir='rl';
i=d-1;
cry_left=cry(1:psy(d)-1);
cry_right=cry(psy(d):psy(d+1)-1);
while ( swp <= nswp && not_converged )
       %As usual, compute (local) B-matrix in the TT-format
       cra1=cra(psa(i):psa(i+1)-1); 
       cra2=cra(psa(i+1):psa(i+2)-1);
       px1=phi{i}; px1=reshape(px1,[ra(i),ry(i),ry(i)]); 
       %px1=permute(px1,[1,2,3]); 
       px1=permute(px1,[1,3,2]);
       px1=reshape(px1,[ra(i),ry(i)*ry(i)]);
       px2=phi{i+2}; 
       px2=reshape(px2,[ra(i+2),ry(i+2),ry(i+2)]);
       %px2=permute(px2,[1,3,2]);
       %Compute the local matrix just by putting px1 into cra1, px2 into
       %cra2
       cra2=reshape(cra2,[ra(i+1)*n(i+1)*n(i+1),ra(i+2)]);
       px2=reshape(px2,[ra(i+2),ry(i+2)*ry(i+2)]);
       b2=cra2*px2; % 
       b1=px1.'*reshape(cra1,[ra(i),numel(cra1)/ra(i)]);
       b1=reshape(b1,[ry(i),ry(i),n(i),n(i),ra(i+1)]); b1_save=b1; %Save for phi computations
       b1=permute(b1,[1,3,2,4,5]); 
       b1=reshape(b1,[ry(i)*n(i),ry(i)*n(i),ra(i+1)]);
       b2=reshape(b2,[ra(i+1),n(i+1),n(i+1),ry(i+2),ry(i+2)]);
       b2_save=b2; %Save for phi computations
       b2=permute(b2,[2,4,3,5,1]); 
       b2=reshape(b2,[n(i+1)*ry(i+2),n(i+1)*ry(i+2),ra(i+1)]);
       mm=cell(2,1); mm{1}=b1; mm{2}=b2;
       cry=[cry_left,cry_right];
       y.r=ry;
       y.core=cry;
       y.ps=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
       %mm1=tt_matrix(mm);
       %Now setup the initial guess: i the core  
       if ( strcmp(dir,'rl') ) %The block index is in the second core
         pos=numel(cry_left);
         cry1=cry_left(pos-ry(i)*n(i)*ry(i+1)+1:pos);
         cry2=cry_right(1:ry(i+1)*n(i+1)*k*ry(i+2));
         cry1=reshape(cry1,[ry(i)*n(i),ry(i+1)]);
         cry2=reshape(cry2,[ry(i+1),n(i+1)*k*ry(i+2)]);
         w=cry1*cry2; w=reshape(w,[ry(i),n(i),n(i+1),k,ry(i+2)]);
         w=permute(w,[1,2,3,5,4]); w=reshape(w,[numel(w)/k,k]);
       else %The block index is in the first core
           pos=numel(cry_left); 
           cry1=cry_left(pos-ry(i)*n(i)*k*ry(i+1)+1:pos);
           cry2=cry_right(1:ry(i+1)*n(i+1)*ry(i+2));
           cry1=reshape(cry1,[ry(i)*n(i)*k,ry(i+1)]);
           cry2=reshape(cry2,[ry(i+1),n(i+1)*ry(i+2)]);
           w=cry1*cry2; w=reshape(w,[ry(i),n(i),k,n(i+1),ry(i+2)]);
           w=permute(w,[1,2,4,5,3]); w=reshape(w,[numel(w)/k,k]);
       end
       %Now run the eigenvalue solver
        bw=bfun(mm,w); ev=bw'*w; 
           er0=norm(bw-w*ev,'fro')/norm(w,'fro');
       if ( size(w,1) >= max(5*k,msize) )
           [wnew,ev,fail_flag]=lobpcg(w,@(x) bfun(mm,x),eps,max_l_steps);
       else
          fm=full(tt_matrix(mm)); 
          [v,dg]=eig(fm);
          ev=diag(dg);
          [ev,ind]=sort(ev,'ascend');
          v=v(:,ind);
          wnew=v(:,1:k);
          ev=ev(1:k);
       end
       er1=norm(bfun(mm,wnew)-wnew*diag(ev),'fro')/norm(wnew,'fro');

       fv=sum(ev); %The functional we minimize;
       fprintf('sweep=%d block=%d fv=%f loc solve=%3.2e old_solve=%3.2e \n',swp,i,fv,er1,er0);
        cry=[cry_left,cry_right];
        y.core=cry;
        y.r=ry; 
        y.d=d;
        
       if ( strcmp(dir,'rl') ) %Implant the auxiliary core into the i-th core
           %(a1,i1,a2,a2,i2*g,a3)-> (a1,i1*g,a2,a2,i2,a3)
           %Delete old block from the core_left, add new block to the core
           %right
           wnew=reshape(wnew,[ry(i),n(i),n(i+1),ry(i+2),k]);
           wnew=permute(wnew,[1,2,5,3,4]); wnew=reshape(wnew,[ry(i)*n(i)*k,n(i+1)*ry(i+2)]);
           [u,s,v]=svd(wnew,'econ'); s=diag(s); 
       
           %Truncation block
           rnew=my_chop2(s,eps*norm(s)); 
           u=u(:,1:rnew); s=s(1:rnew); v=v(:,1:rnew);% v=v';
           u=u*diag(s); %u has to be reshaped 

           %Random restart block
           radd=min(kick_rank,size(v,1)-rnew);
           rnew=rnew+radd;
           if ( radd >  0 )
             vr=randn(size(v,1),radd);
             ur=zeros(size(u,1),radd); 
             %Orthogonalize vr to v by Golub-Kahan reorth
             mvr=v'*vr; vnew=vr-v*mvr; 
             reort_flag=false;
             for j=1:radd
                if ( norm(vnew(:,j)) <= 0.5*norm(vr(:,j)))
                   reort_flag=true;
                end
             end
             if (reort_flag)
                 sv=v'*vnew;
                 %mvr=mvr+v'*vnew; 
                 vnew=vnew-v*sv; 
             end
             [vnew,~]=qr(vnew,0); 
             
             v=[v,vnew];
             u=[u,ur];
             %norm(u*v'-u1*v1');
             %keyboard;
             %keyboard;
           end
           v=v';
           
           
           
           
           %Memory stuff
           pos=numel(cry_left);
           cry_left(pos-ry(i)*n(i)*ry(i+1)+1:pos)=[];
           cry_right(1:ry(i+1)*n(i+1)*k*ry(i+2))=[]; %Delete the top core from cry_right
           cry_right=[u(:)', v(:)',cry_right]; %Simply add new block to cry_right
           ry(i+1)=rnew;
           
           %Recalculate phi block; we need to recalculate phi{i+1} here
           %using phi{i+2} 
           %px2(ra(i+2),ry(i+2),ry(i+2))*cra2(ra(i+1),n(i),m(i),ra(i+2)*v(ry(i+1),n(i+1),ry(i+2))*v(ry(i+1),n(i+1),ry(i+2))
           %we already have b2 ---  (b2_save)
           %(ra(i+1),n(i+1),n(i+1),ry(i+2),ry(i+2)), convolve over the 
           %n(i+1),ry(i+2)
           phx=reshape(b2_save,[ra(i+1),n(i+1),n(i+1),ry(i+2),ry(i+2)]); 
           phx=permute(phx,[2,4,1,3,5]);
           %phx=permute(phx,[3,4,1,2,5]);
           %phx=permute(phx,[3,5,1,2,4]);
           % phx=permute(phx,[2,5,1,3,4]);
           phx=reshape(phx,[n(i+1)*ry(i+2),ra(i+1)*n(i+1)*ry(i+2)]);
           v0=v;
           v0=reshape(v0,[ry(i+1),n(i+1)*ry(i+2)]);
           phx=v0*phx;
           phx=reshape(phx,[ry(i+1),ra(i+1),n(i+1),ry(i+2)]);
           phx=permute(phx,[3,4,1,2]); 
           phx=reshape(phx,[n(i+1)*ry(i+2),ry(i+1)*ra(i+1)]);
           phx=v0*phx;
           phx=reshape(phx,[ry(i+1),ry(i+1),ra(i+1)]);
           %phx=permute(phx,[3,2,1]); 
           phx=permute(phx,[3,2,1]);
           phi{i+1}=phx; 
       else %Implant the auxiliary core from the left block to the right block
  
           wnew=reshape(wnew,[ry(i),n(i),n(i+1),ry(i+2),k]); 
           wnew=permute(wnew,[1,2,3,5,4]); wnew=reshape(wnew,[ry(i)*n(i),n(i+1)*k*ry(i+2)]);
           [u,s,v]=svd(wnew,'econ'); s=diag(s);
           
           %Truncation block
           rnew=my_chop2(s,eps*norm(s));
           u=u(:,1:rnew); s=s(1:rnew); v=v(:,1:rnew); v=v*diag(s); 
           
           %Random restart block
           radd=min(kick_rank,size(u,1)-rnew);
           rnew=rnew+radd;
           if ( radd >  0 )
             vr=zeros(size(v,1),radd);
             ur=randn(size(u,1),radd); 
             %Orthogonalize vr to v by Golub-Kahan reorth
             mur=u'*ur; unew=ur-u*mur; 
             reort_flag=false;
             for j=1:radd
                if ( norm(unew(:,j)) <= 0.5*norm(ur(:,j)))
                   reort_flag=true;
                end
             end
             if (reort_flag)
                 sv=u'*unew;
                 %mvr=mvr+v'*vnew; 
                 unew=unew-u*sv; 
             end
             [unew,~]=qr(unew,0); 
             
             u=[u,unew];
             v=[v,vr];
           end
           v=v';
           
           
           
           cry_right(1:ry(i+1)*n(i+1)*ry(i+2))=[];
           pos=numel(cry_left);
           cry_left(pos-ry(i)*n(i)*k*ry(i+1)+1:pos)=[];
           cry_left=[cry_left,u(:)',v(:)'];
           
           
           
           
           
           
           
           ry(i+1)=rnew;
           
           
           
           
           %Recalculate phi block; we need to recalculate phi{i+1} using
           %phi{i}
           phx=reshape(b1_save,[ry(i),ry(i),n(i),n(i),ra(i+1)]);
           u0=u; u0=reshape(u0,[ry(i)*n(i),ry(i+1)]);
           phx=permute(phx,[1,3,5,2,4]); 
           %phx=permute(phx,[2,4,5,1,3]); 
           %phx=permute(phx,[1,4,5,2,3]);
           %phx=permute(phx,[2,3,5,1,4]);
           phx=reshape(phx,[ry(i)*n(i)*ra(i+1),ry(i)*n(i)]);
           phx=phx*u0; 
           phx=reshape(phx,[ry(i),n(i),ra(i+1),ry(i+1)]);
           phx=permute(phx,[3,4,1,2]); 
           phx=reshape(phx,[ra(i+1)*ry(i+1),ry(i)*n(i)]);
           phx=phx*u0; phx=reshape(phx,[ra(i+1),ry(i+1),ry(i+1)]);
           phx=permute(phx,[1,2,3]);
           phi{i+1}=phx;
       end
   %Choose the next direction block; now implement the simple case
   if ( strcmp(dir,'rl') )
     if ( i > 1 )
       i=i-1;
     else %Change direction %The last optimization was for (1,2) core 
       dir='lr';
       %One block should go from cry_right to cry_left
       cry_left=cry_right(1:ry(1)*n(1)*ry(2)*k); %This seems correct
       cry_right(1:ry(1)*n(1)*ry(2)*k)=[];
       swp=swp+1;
     end
   else
    if ( i < d-1 )
        i=i+1;
     else
       dir='rl';
       pos=numel(cry_left);
       cry_right=cry_left(pos-ry(d)*n(d)*ry(d+1)*k+1:pos); 
       cry_left(pos-ry(d)*n(d)*ry(d+1)*k+1:pos)=[];
       swp=swp+1;
       %One block should go from cry_left to cry_right (?) --- seems no :)
     end
   end
end

    function [x]=bfun(a,x)
    %[Y]=BFUN(A,X,K)
    %a is given as U(i1,j1,s)*V(i2,j2,s)
    %\sum_{j1,j2,s} U(i1,j1,s)*V(i2,j2,s)*X(j1,j2,q)
    re=size(x,2);
    ul=a{1}; vl=a{2};
    n1=size(ul,1);m1=size(ul,2); ral=size(ul,3);
    n2=size(vl,1); m2=size(vl,2);
    ul=permute(ul,[1,3,2]); ul=reshape(ul,[numel(ul)/m1,m1]);
    x=reshape(x,[m1,numel(x)/m1]);
    x=ul*x; %x is i1xsxj2xq s,j2
    x=reshape(x,[n1,ral,m2,re]);
    x=permute(x,[3,2,1,4]); 
    x=reshape(x,[m2*ral,n1*re]);
    vl=reshape(vl,[n2,m2*ral]); 
    x=vl*x; %is n2*n1*k
    x=reshape(x,[n2,n1,re]); 
    x=permute(x,[2,1,3]); 
    x=reshape(x,[numel(x)/re,re]);
    return
    end

return
end