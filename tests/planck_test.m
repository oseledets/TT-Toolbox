dpx = 8; % phys. dims for x
nn=4;
ddd=5;
zzz=0.1;
xx=cell(nn,3); %The space is three-dimensional
a=-10;
b=10;
np=2^(dpx);
h=(b-a)/(np+1);
e=tt_ones(dpx*nn,2); e=tt_tensor(e);
for i=1:nn
  xx0=tt_xtr(dpx,nn,i);
  xx0=a*e+(xx0+e)*h;
  xx0=round(xx0,1e-8);
  xx{i,1}=kron(kron(xx0,e),e);
  
  xx{i,2}=kron(kron(e,xx0),e);
  
  xx{i,3}=kron(kron(e,e),xx0);
end

%Now be careful and generate the TRUE potential
z=diag(ones(nn-1,1),1);
e0=eye(nn);
lp=2*e0-z-z';
[vec,ev]=eig(lp);
ev=diag(ev);
[ev,ind]=sort(ev,'descend');
vec=vec(:,ind);
V=vec';%;*diag(sqrt(ev));
%inv(V)*lp*inv(V')
%return
%norm(lp-vec*diag(ev)*vec')
%return
%V=zeros(nn);
%for i=1:nn
%  for j=1:nn
%    V(i,j) = sin(i*j*pi/(nn+1));
%  end
%end

%Create the transformed coordinates;
xxT=cell(nn,3);
for i=1:nn
  xxT{i,1}=xx{1,1}*V(i,1);
  xxT{i,2}=xx{1,2}*V(i,1);
  xxT{i,3}=xx{1,3}*V(i,1);
  
  for j=2:nn
    xxT{i,1}=xxT{i,1}+xx{j,1}*V(i,j);
    xxT{i,2}=xxT{i,2}+xx{j,2}*V(i,j);
    xxT{i,3}=xxT{i,3}+xx{j,3}*V(i,j);
  end
  xxT{i,1}=round(xxT{i,1},1e-8);
  
  xxT{i,2}=round(xxT{i,2},1e-8);
  
  
  xxT{i,3}=round(xxT{i,3},1e-8);
end
%Compute pairwise distances
ruv=cell(nn,nn);
for i=1:nn
  for j=i:nn
%      pp1=xxT{i,1};
%      pp2=xxT{i,2};
%      pp3=xxT{i,3};
%      
%      for k=i+1:j
%         pp1=pp1+xxT{k,1};
%         pp1=round(pp1,1e-10);
%         pp2=pp2+xxT{k,2};
%         pp2=round(pp2,1e-10);
%         pp3=pp3+xxT{k,3};
%         pp3=round(pp3,1e-10);
%      end
     pp1=xxT{i,1}-xxT{j,1};
     
     pp2=xxT{i,2}-xxT{j,2};
     
     pp3=xxT{i,3}-xxT{j,3};
     
     ruv{i,j}=pp1.^2+pp2.^2+pp3.^2;
     ruv{i,j}=round(ruv{i,j},1e-8);
     ruv{j,i}=ruv{i,j};
  end
end
eexp=cell(nn,nn);
for i=1:nn
    for j=i:nn
       eexp{i,j}=funcross(ruv{i,j},@(x) exp(-0.5*x/ddd.^2),1e-7,ruv{i,j},20);
       eexp{i,j}=round(eexp{i,j},1e-7);
       eexp{j,i}=eexp{i,j};
    end
end
%We have to sum them all
pot=1e-16*eexp{1,1}/norm(eexp{1,1}); pot=round(pot,2);
for i=1:nn
  for j=1:nn
      if ( i ~= j )
      pot=pot+eexp{i,j};
      pot=round(pot,1e-7);
      end
  end
end
%And the harmonic part
% pot=pot*zzz*0.5;
% for i=1:nn
%    pot=pot+xxT{i}.^2;
%    pot=round(pot,1e-7);
% end
%return
 %mm1=funcross(rr,@(x) exp(-0.5*x/ddd.^2),2*h^2,mm,100);
 %if ( norm(mm-mm1)/norm(mm) > h )
 %  fprintf('Breakdown, sir! \n');
 %end
%rr=cell(nn,nn);
%rr{1}=xx{1};
%for i=2:n
% rr{i}=rr{i-1}+xx{i};
% rr{i}=round(rr{i},2);
%end
rr=xx{1}-xx{2}; rr=rr.^2;
rr=round(rr,1e-8);
mm=funcross(rr,@(x) exp(-0.5*x/ddd.^2),h^2,rr,100);
 mm1=funcross(rr,@(x) exp(-0.5*x/ddd.^2),2*h^2,mm,100);
 if ( norm(mm-mm1)/norm(mm) > h )
   fprintf('Breakdown, sir! \n');
 end

% rr=xx{1};
% for i=2:nn
%   rr=rr+xx{i};
%   rr=round(rr,1e-12);
% end
% %rr=round(rr,1e-8);
% sgn1=funcross(rr,@(x) sign(x),1e-9,rr,100);
% sgn2=funcross(rr,@(x) sign(x),1e-8,rr,100);
% if ( norm(sgn1-sgn2)/norm(sgn2) > 1e-4 )
%   fprintf('Breakdown, sir! \n');
% end
% %sgn1=full(sgn1,2*ones(1,2*dpx));
% %prm=1:2*dpx; prm=reshape(prm,[dpx,2]); prm=prm';
% %prm=reshape(prm,[1,2*dpx]);
% %sgn1=ipermute(sgn1,prm);
% %sgn1=reshape(sgn1,np,np);
