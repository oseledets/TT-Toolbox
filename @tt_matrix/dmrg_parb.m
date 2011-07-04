function [y]=dmrg_parb(mat,a,f,eps,y0,rmax,nswp,verb)
%[y]=dmrg_parb(mat,a,rhs,eps,[y0],[rmax],[nswp],[ver])
%This is not an east stuff.
%It solve parameter-dependent problem A(y)u(y) = rhs(y)
%With accuracy eps.
%A(y) is given as a sum mat{s}*a(s,P), and a is stored
%as a row TT-tensor (i.e. r(1) = r
%mat is stored as 1D TT-matrix with r(2) = r

%The storage is a delicate issue. We will stick to the
%separate storage of the "Physical" part and of the diagonal part
%We will focus on the "maxvol" part.


%Parameters section
kick_rank=5;

%This is the "physical dimension" size

n=a.n;
d=a.d; %And this is for the parametric part
%We seek for the solution as JxY1xQ1xQ1xY2xQ2
k=rank(mat.tt,2); %This is "the block size"
N=size(mat,1);

matP=full(mat.tt); 
matP=reshape(matP,[N*N,k]);
%matP=full(mat);
n=[N;n];
if ( nargin < 4 || isempty(y0) )
   y=tt_tensor;
   y.d=d+1;
   ry=k*ones(1,d); ry(1)=k; ry(d+1)=1; ry=ry';
   y.r=ry;
   y.n=n;
   y.ps=cumsum([1;y.n.*y.r(1:y.d).*y.r(2:y.d+1)]);
   psy=y.ps;
   sz=psy(d+1)-1;
   cr=randn(1,sz);
   %cr=ones(sz,1);
   y.core=cr;
else 
   y=y0;
end
if ( nargin < 5 || isempty(rmax) )
  rmax=1000;
end
if ( nargin < 6 || isempty(nswp) )
  nswp=10;
end
if ( nargin < 7 || isempty(verb) )
  verb=true;
end
%Matrix details
cra=a.core;
ra=a.r;
psa=a.ps;
ry=y.r;
psy=y.ps;
cry=y.core;
crf=f.core;
psf=f.ps;
fP=crf(psf(1):psf(2)-1); fP=reshape(fP,[N,numel(fP)/N]);
kf=size(fP,2);
crf(psf(1):psf(2)-1)=[];
%The user will never know what happens here
rf=f.r;
rf=rf(2:d+2);
nf=f.n;
nf=nf(2:d+1);
psf=cumsum([1;nf.*rf(1:d).*rf(2:d+1)]);


%Right hand side details --- we would like to split 
%the physical part and the parametric part here also,
%but we do not want to bother the reader with such amount
%of useless parameters, so we will do it ourselves(?)
%NO! We can just compute phi up to i=2 and that is all




%The situation is "kinda" similar; however, we start from right-to-left 
%qr & maxvol of Y; The "matrix" core is like ra(i)xn(i)xm(i)*ra(i+1)
%That means, that the "left" matrix is ry(i)xkxk, and the right is just
%ra(i+1)
%Actually, we need for a specific index subset (y1,...,yk,J,yk+1...yP)
%Compute the matrix; it is more (or less) reduces to the computation
%of the coefficient; i.e. B(i,j,Q)*Phi(Q,ALLINDICES) -> B(I,J,ALLINDICES)
%B(i,j,Q1)*C(Q1,y1,Q2)*C(Q2,y2,Q3)*C(Q3,y3,Q4)*C(Q4,y4) 
%B(i,j,Q1)*P(Q1,Q3,RY1)*C(Q3,y3,Q4)*P(Q4,RY2,Q5)
%B(i,j,Q3,RY1)
yold=y;
rm=1;
pos1=psy(d+2)-1;
phi=cell(d+1,1);
phi{d+1}=1;
phif=cell(d+1,1);
phif{d+1}=1;

%Right-to-left qr & maxvol
for i=d+1:-1:2
    %fprintf('i=%d psy(%d)=%d psy(%d)=%d ry(%d)=%d n(%d)=%d ry(%d)=%d \n',i,i,psy(i),i+1,psy(i+1),i,ry(i),i,n(i),i+1,ry(i+1));
    cr=cry(psy(i):psy(i+1)-1);
    cr=reshape(cr,[ry(i)*n(i),ry(i+1)]);
    cr=cr*rm; 
    %cry1=cry;
    %cry1(psy(i):psy(i+1)-1)=cr(:);
    %y.core=cry1;
    %norm(nm*y-f)
    %keyboard

    ry(i+1)=size(cr,2); 
    cr=reshape(cr,[ry(i),n(i)*ry(i+1)]);
    cr=cr.'; %It seems to be needed
    [u,rm]=qr(cr,0); 
    indu=maxvol2(u); 
    r1=u(indu,:);
    u=u/r1;
    rm=r1*rm; rm=rm.';
    rnew=size(u,2);
    u=u.'; %n(i)*ry(i+1)xry(i)->ry(i)xn(i)xry(i+1) 
    %ry(i+1)=rnew;
    
    cry(pos1-ry(i+1)*n(i)*rnew+1:pos1)=u(:); 
    pos1=pos1-ry(i+1)*n(i)*rnew; %Mda
    %With this new core compute new phi matrix; we better need raxry
    
    crm=cra(psa(i-1):psa(i)-1);
    phx=phi{i};
    phx=reshape(phx,[ra(i),ry(i+1)]);
    crm=reshape(crm,[ra(i-1)*n(i),ra(i)]);
    crm=crm*phx; %crm is ra(i)xn(i)xry(i+1); u is rnewxn(i)xry(i+1);
    crm=reshape(crm,[ra(i-1),n(i)*ry(i+1)]); 
    phi{i-1}=crm(:,indu);
    crh=crf(psf(i-1):psf(i)-1); 
    crh=reshape(crh,[rf(i-1)*n(i),rf(i)]);
    phf=reshape(phif{i},[rf(i),ry(i+1)]);
    crh=crh*phf;
    crh=reshape(crh,[rf(i),n(i)*ry(i+1)]);
    phif{i-1}=crh(:,indu);
end
cr=cry(psy(1):psy(2)-1);
cr=reshape(cr,[ry(1)*n(1),ry(2)]);
cr=cr*rm; 
ry(2)=size(cr,2); 
cry(pos1-ry(1)*n(1)*ry(2)+1:pos1)=cr(:);
pos1=pos1-ry(1)*n(1)*ry(2)+1; %Fix it to the start
%And now truncate cry
cry=cry(pos1:end);

%Test that here all is okey.
psy=cumsum([1;(n).*ry(1:y.d).*ry(2:y.d+1)]);
y.core=cry;
y.r=ry;
y.ps=psy;
%norm(nm*y-f)
%keyboard


%Now we have to implant the "physical dimension" into the parametric 
%part for further floating
psy=cumsum([1;(n).*ry(1:d+1).*ry(2:d+2)]);
cr1=cry(psy(1):psy(2)-1);
cr2=cry(psy(2):psy(3)-1); 
cr1=reshape(cr1,[numel(cr1)/ry(2),ry(2)]);
cr2=reshape(cr2,[ry(2),numel(cr2)/ry(2)]);
cr_new=cr1*cr2; %crnew is now NxPxRY(3)
cr_new=reshape(cr_new,[N,n(2),ry(3)]);
%cr_new=permute(cr_new,[2,1,3]);
cry(1:psy(3)-1)=[];
cry=[cr_new(:);cry]; %Well now it is prepared
ry=[N;ry(3:end)]; %Now the physical dimension is a hanging rank
n=n(2:end); 
psy=cumsum([1;(n).*ry(1:d).*ry(2:d+1)]);
y.core=cry;
y.r=ry;
y.n=n;
y.d=d;
y.ps=psy;
%Start the main iteration (borrow the idea from dmrg_eigb code)


psy=y.ps;
cry=y.core;
swp=1;
not_converged=true;
dir='lr';
i=1;
cry_left=cry(1:psy(2)-1); %Bullshit
cry_right=cry(psy(2):psy(d+1)-1);
%P(1)xRA(2) select "good rows" out of it
%First PHI matrix should be of size QxRY(1)xRA(1); what if RY(1) = 1?
%We know that the considered core is Q times larger
ry(1)=1; %This should be done
%We computed phi up to i=1 the last one was the first parametric
%mode i.e. (NXN x Q) x [Q x RY] and that is the true thing;
%We then can turn it to (NxN)xRY and start with K=RY
%B(I,J,Q)xDELTA(Q,Q')*[(Q'xP1xRA2)x(RA2xP2x)],S)-> 
%B(I,J,Q)

%Stupidity: gather the local matrix
phi{1}=eye(k); %Heh, this seems ok now :) Maybe the bug is 
%phi{1}=diag(phi{1}(:));
%here --- we have had phi{1} of size ra(1)xry(1); 
phif{1}=eye(kf); %The same
%phif{1}=diag(phif{1}(:));
%B(I,J,Q)*Z(Q)*CA(Q,P1,RA2)*CA(RA2,P2,RA3)*CA(RA3,S)

while ( swp <= nswp && not_converged )
   %Gather the local matrix.( The left phi matrix always includes
   %the  "matrix" rank)
   %Let us look at the summation
   %it is B(i,j,Q)*PHI(Q,RY1,N1,N2,RY2) HEY! What happens for the first
   %core?
   %Now how the matrix PHI is gathered
   %PHI(Q,RY1,RA1)*B1(RA1,N1,RA2)*B2(RA2,N2,RA3)*PHI(RA3,RY3)
   %Seems rather simple...not yet :(((
   ph1=phi{i};
   ph2=phi{i+2}; 
   cra1=cra(psa(i):psa(i+1)-1);
   cra2=cra(psa(i+1):psa(i+2)-1);
   ph1=reshape(ph1,[k*ry(i),ra(i)]); %
   cra1=reshape(cra1,[ra(i),n(i)*ra(i+1)]);
   ph1=ph1*cra1; %ph1 is k*ry(i)*n(i)*ra(i+1);
   b1=ph1; %Save for phi calculations
   b1=reshape(b1,[numel(b1)/ra(i+1),ra(i+1)]);
   cra2=reshape(cra2,[ra(i+1)*n(i+1),ra(i+2)]);
   ph2=reshape(ph2,[ra(i+2),ry(i+2)]);
   b2=cra2*ph2;
   b2=reshape(b2,[ra(i+1),numel(b2)/ra(i+1)]);
   ph1=b1*b2;
   ph1=reshape(ph1,[k,numel(ph1)/k]);
   B=matP*ph1; %B is NxNxry(i)*n(i)*n(i+1)*ry(i+2)
   py=ry(i)*n(i)*n(i+1)*ry(i+2);
   B=reshape(B,N,N,py);
   %Now form the rhs. It it is somewhat more complicated. 
   %For simplicity, we will split the stuff as for the matrix 
   %in the same fashion
   ph1=phif{i};
   ph2=phif{i+2}; 
   crf1=crf(psf(i):psf(i+1)-1);
   crf2=crf(psf(i+1):psf(i+2)-1);
   ph1=reshape(ph1,[kf*ry(i),rf(i)]); %
   crf1=reshape(crf1,[rf(i),n(i)*rf(i+1)]);
   ph1=ph1*crf1;
   bf1=ph1; %Save for phi calculations
   bf1=reshape(bf1,[numel(bf1)/rf(i+1),rf(i+1)]);
   crf2=reshape(crf2,[rf(i+1)*n(i+1),rf(i+2)]);
   ph2=reshape(ph2,[rf(i+2),ry(i+2)]);
   bf2=crf2*ph2;
   bf2=reshape(bf2,[rf(i+1),numel(bf2)/rf(i+1)]);
   ph1=bf1*bf2;
   ph1=reshape(ph1,[kf,numel(ph1)/kf]);
   fB=fP*ph1; %fB is Nxry(i)*n(i)*n(i+1)*ry(i+2)
   py=ry(i)*n(i)*n(i+1)*ry(i+2);
   fB=reshape(fB,N,py);
   solB=zeros(N,py);
   %tic;
   %for j=1:py
   %  solB(:,j)=reshape(B(:,:,j),N,N) \ fB(:,j); %Bug fixed
   %end
   %t1=toc;
   %tic;
   solB=parmldivide(B,fB);
   %t2=toc;
   %fprintf('save=%3.1f \n',t1/t2);
   %toc;
   %if ( py > 1000 )
   %  keyboard
   %end
   %solB is Nxry(i)*n(i)*n(i+1)*ry(i+2)
   solB=reshape(solB,[N,ry(i),n(i),n(i+1),ry(i+2)]);
   %Compute the previous solution
   if ( strcmp(dir,'lr') )
     pos=numel(cry_left);
     w1=cry_left(pos-ry(i)*N*n(i)*ry(i+1)+1:pos);
     w2=cry_right(1:ry(i+1)*n(i+1)*ry(i+2));
     w1=reshape(w1,[numel(w1)/ry(i+1),ry(i+1)]);
     w2=reshape(w2,[ry(i+1),numel(w2)/ry(i+1)]);
     w=w1*w2; w=reshape(w,[ry(i),N,n(i),n(i+1),ry(i+2)]);
     w=permute(w,[2,1,3,4,5]); 
   elseif (strcmp(dir,'rl') )
          pos=numel(cry_left);
     w1=cry_left(pos-ry(i)*n(i)*ry(i+1)+1:pos);
     w2=cry_right(1:ry(i+1)*N*n(i+1)*ry(i+2));
     w1=reshape(w1,[numel(w1)/ry(i+1),ry(i+1)]);  %Final step is like [ry(i)xNxn(i)]n(i+1)xry(i+2) -> ry(i)xn(i)x [Nxn(i+1)*ry(i+2)]
     w2=reshape(w2,[ry(i+1),numel(w2)/ry(i+1)]);
     w=w1*w2; 
     %w=reshape(w,[ry(i),n(i),N,n(i+1),ry(i+2)]);
     w=reshape(w,[ry(i),n(i),n(i+1),N,ry(i+2)]);
     w=permute(w,[4,1,2,3,5]); 
     
     %w=permute(w,[3,1,2,4,5]); 
   end
   %if ( i == 7 ) 
   %  keyboard;
   %end
   er1=norm(w(:)-solB(:))/norm(w(:));
   if ( verb )
  fprintf('sweep=%d block=%d  error=%3.2e \n',swp,i,er1);
   end
   %And now add the splitting of solB part + recomputation of phi & phif
   %Memory stuff
   if ( strcmp(dir,'lr') ) %Implant the auxiliary core from the left block to the right block 
       %solB is Nxry(i)*n(i)xn(i+1)*ry(i+2) %ry(i+1)*N*n(i+1)*ry(i+2)
       solB=permute(solB,[2,3,1,4,5]); solB=reshape(solB,[ry(i)*n(i),N*n(i+1)*ry(i+2)]);
       [u,s,v]=svd(solB,'econ');
       s=diag(s);
       rnew=my_chop2(s,eps*norm(s));
       u=u(:,1:rnew); s=s(1:rnew); v=v(:,1:rnew); v=v*diag(s);      
       ur=randn(size(u,1),kick_rank);
         %Orthogonalize ur to u by Golub-Kahan reorth
         u=reort(u,ur);
         radd=size(u,2)-rnew; 
         if ( radd > 0 )
           vr=zeros(size(v,1),radd);
           v=[v,vr];
         end
         rnew=rnew+radd;
         indu=maxvol2(u); 
         r1=u(indu,:);
         u=u/r1;
         v=v*r1';
         v=v';
         %Memory stuff
          
         cry_right(1:ry(i+1)*n(i+1)*ry(i+2))=[];
         pos=numel(cry_left);
         cry_left(pos-ry(i)*n(i)*N*ry(i+1)+1:pos)=[];
         cry_left=[cry_left,u(:)',v(:)'];
         ry(i+1)=rnew;
           
           
           
         %Recalculate phi block; we need to recalculate phi{i+1} using
         %phi{i}
         %b1 is Qxry(i)*n(i)*ra(i+1);
         b1=reshape(b1,k,ry(i)*n(i),ra(i+1));
         phi{i+1}=b1(:,indu,:);
         bf1=reshape(bf1,kf,ry(i)*n(i),rf(i+1));
         phif{i+1}=bf1(:,indu,:);
         %That is all, folks!
       
    elseif ( strcmp(dir,'rl') ) %Implant the auxiliary core into the i-th core
      %solB is [Nxry(i)xn(i)xn(i+1)xry(i+2)] it was like
      %ry(i)*n(i)*n(i+1)*N*ry(i+2) for implanting 
       solB=permute(solB,[2,3,1,4,5]); solB=reshape(solB,[ry(i)*n(i)*N,n(i+1)*ry(i+2)]);
      %Truncation block
      [u,s,v]=svd(solB,'econ');
      s=diag(s); rnew=my_chop2(s,norm(rnew)*eps);
      
      u=u(:,1:rnew); s=s(1:rnew); v=v(:,1:rnew);% v=v';
      u=u*diag(s); %u has to be reshaped 
       
      vr=randn(size(v,1),kick_rank);   
      v=reort(v,vr);
      radd=size(v,2)-rnew; 
      if ( radd > 0 )
         ur=zeros(size(u,1),radd);
         u=[u,ur];
      end
      rnew=rnew+radd;
      indv=maxvol2(v); 
      r1=v(indv,:);
      v=v/r1;
      u=u*r1';
      v=v';
      %Memory stuff
      pos=numel(cry_left);
      cry_left(pos-ry(i)*n(i)*ry(i+1)+1:pos)=[];
      cry_right(1:ry(i+1)*n(i+1)*N*ry(i+2))=[]; %Delete the top core from cry_right
      cry_right=[u(:)', v(:)',cry_right]; %Simply add new block to cry_right
      ry(i+1)=rnew;
      b2=reshape(b2,[ra(i+1),n(i+1)*ry(i+2)]);
      b2=b2(:,indv);
      phi{i+1}=b2;
      bf2=reshape(bf2,[rf(i+1),n(i+1)*ry(i+2)]);
      bf2=bf2(:,indv);
      phif{i+1}=bf2;
      
      %New block 
      
   end
   %Now the test what to do in the end
      %Choose the next direction block; now implement the simple case
   if ( strcmp(dir,'rl') )
     if ( i > 1 )
       i=i-1;
     else %Change direction %The last optimization was for (1,2) core 
       dir='lr';
       %One block should go from cry_right to cry_left
       cry_left=cry_right(1:ry(1)*n(1)*ry(2)*N); %This seems correct
       cry_left=reshape(cry_left,[ry(1),n(1),N,ry(2)]);
       cry_left=permute(cry_left,[1,3,2,4]);
       cry_left=cry_left(:);
       cry_right(1:ry(1)*n(1)*ry(2)*N)=[];
       swp=swp+1;
       
     end
   else
    if ( i < d-1 )
        i=i+1;
     else
       dir='rl';
       pos=numel(cry_left);
       cry_right=cry_left(pos-ry(d)*n(d)*ry(d+1)*N+1:pos); 
       cry_right=reshape(cry_right,[ry(d),N,n(d),ry(d+1)]);
       cry_right=permute(cry_right,[1,3,2,4]);
       cry_right=cry_right(:);
       cry_left(pos-ry(d)*n(d)*ry(d+1)*N+1:pos)=[];
       %One block should go from cry_left to cry_right (?) --- seems no :)
     end
   end
end
%All we have to do is to join cry_left and cry_right
cry=[cry_left(:);cry_right(:)];
%We have to reshape it back to the old structure; 
ry(1)=N; 
psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
cr1=cry(psy(1):psy(2)-1);
cr1=reshape(cr1,[N,n(1),ry(2)]);
%cr1=permute(cr1,[2,1,3]);
%cr1=reshape(cr1,[N,n(1),ry(2)]);
cr1=reshape(cr1,[N,numel(cr1)/N]);
[u,s,v]=svd(cr1,'econ');
cry(psy(1):psy(2)-1)=[];
u=u*s; v=v'; rfin=size(u,2);
cry=[u(:);v(:);cry(:)];
ry(1)=[];
ry=[1;rfin;ry(:)];
n=[N;n];
d=d+1;
psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
y.ps=psy;
y.r=ry;
y.n=n;
y.d=d;
y.core=cry;
return
end
