function [tt]=get_projector(d,a,b)
%[TT]=GET_PROJECTOR(D,A,B)
%Simple 1D linear interpolator from coarse (XX1) to fine (XX2) 
%uniform grids on [A,B]
%grid with 2^(D+1) nodes
%in QTT format
n1=2^d;
n2=2^(d+1);

h1=(b-a)/(n1+1); h2=(b-a)/(n2+1);
xx1=a+(1:n1)*h1; xx2=a+(1:n2)*h2;
prj=zeros(n2,n1); 
for i=1:n1-1
   prj(2*i,i) = -(xx2(2*i)-xx1(i+1))/(xx1(i+1)-xx1(i));   
   prj(2*i,i+1)=(xx2(2*i)-xx1(i))/(xx1(i+1)-xx1(i));
   prj(2*i+1,i)=-(xx2(2*i+1)-xx1(i+1))/(xx1(i+1)-xx1(i));
   prj(2*i+1,i+1)=(xx2(2*i+1)-xx1(i))/(xx1(i+1)-xx1(i));
end
prj(1,1)=1;
prj(n2,n1)=1;
sz=[2*ones(1,d+1),2*ones(1,d),1];
 prj1=reshape(prj,sz); prm=1:2*(d+1); prm=reshape(prm,[d+1,2]);
prm=prm'; prm=reshape(prm,[1,2*(d+1)]); 
sz2=[4*ones(1,d),2];
prj1=permute(prj1,prm); prj1=reshape(prj1,sz2);
prj_tt=full_to_tt(prj1,1e-8);
tt=tt_vec_to_mat(prj_tt,2*ones(1,d+1),[2*ones(1,d),1]);

return
end