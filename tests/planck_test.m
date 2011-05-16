dpx = 10; % phys. dims for x
nn=2;

xx=cell(nn,1);
a=-10;
b=10;
np=2^(dpx);
h=(b-a)/(np+1);
e=tt_ones(dpx*nn,2); e=tt_tensor(e);
for i=1:nn
  xx{i}=tt_xtr(dpx,nn,i);
  xx{i}=a*e+(xx{i}+e)*h;
  xx{i}=round(xx{i},1e-8);
end

%Now be careful and generate the TRUE potential
V=zeros(nn);
for i=1:nn
  for j=1:nn
    V(i,j) = sqrt(2/nn+1)*sin(i*j*pi/(nn+1));
  end
end
V=V/sqrt(nn+1);
%Create the transformed coordinates;
xxT=cell(nn,1);
for i=1:nn
  xxT{i}=xx{1}*V(i,1);
  for j=2:nn
    xxT{i}=xxT{i}+xx{j}*V(i,j);
  end
  xxT{i}=round(xxT{i},1e-8);
end
w1=xxT{1}.^2; w2=xxT{2}.^2; w12=(xxT{1}+xxT{2}).^2;
ww1=funcross(w1,@(x) exp(-0.5*x/ddd.^2),1e-7,w1,100);

ww12=funcross(w12,@(x) exp(-0.5*x/ddd.^2),1e-7,w12,100);

return
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
