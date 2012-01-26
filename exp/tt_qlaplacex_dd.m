function [tt]=tt_qlaplacex_dd(d)

% returns a rank-5 QTT decomposition of
% Delta_{1}^{-1} \ otimes \ldots \otimes Delta_{D}^{-1},
% Delta_{k} being a discretization of Laplace operator on 2^d(k) points
% uniform grid,
% Dirichlet-Dirichlet boundary conditions being imposed
%
% D=size(d,2) must be >= 1
%
% January 17, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%
% Look for details in the Preprint No. 75, 2010 of
% Max-Planck Institute for Mathematics in the Sciences
% Vladimir A. Kazeev and Boris N. Khoromskij
% On explicit QTT representation of Laplace operator and its inverse
% http://www.mis.mpg.de/publications/preprints/2010/prepr2010-75.html

d=fliplr(d);
D=size(d,2);
tt=cell(sum(d),1);
I=eye(2);
P=[0,1;1,0];
E=ones(2);
F=[1,-1;-1,1];
K=[-1,0;0,1];
L=[0,-1;1,0];


key=0;
for k=1 : D
	for kappa=1 : d(k)
		key=key+1;
		xi=(2^(kappa-1)+1)/(2^kappa+1);
		eta=2^(kappa-2)/(2^kappa+1);
		zeta=xi*(2^(kappa-1)+1)/2^(kappa-1);
		if (kappa == 1)
			tt{key}=zeros(2,2,5);
			tt{key}(:,:,1)=(E+I)/3;
			tt{key}(:,:,2)=2*E;
			tt{key}(:,:,3)=F/18;
			tt{key}(:,:,4)=2*K/3;
			tt{key}(:,:,5)=-2*L/3;
			if (k ~= 1)
				tt{key}=permute(tt{key},[1,2,4,3]);
			end
		elseif (kappa == d(k))
			tt{key}=zeros(2,2,5);
			tt{key}(:,:,1)=I;
			tt{key}(:,:,2)=xi*I/4+zeta*P/4;
			tt{key}(:,:,3)=xi*I-zeta*P;
			tt{key}(:,:,4)=-xi*K/2;
			tt{key}(:,:,5)=zeta*L/2;
		else
			tt{key}=zeros(2,2,5,5);
			tt{key}(:,:,1,1)=I;
			tt{key}(:,:,2,1)=xi*I/4+zeta*P/4;
			tt{key}(:,:,3,1)=xi*I-zeta*P;
			tt{key}(:,:,4,1)=-xi*K/2;
			tt{key}(:,:,5,1)=zeta*L/2;
			tt{key}(:,:,2,2)=2*E;
			tt{key}(:,:,2,3)=2*eta^2*F;
			tt{key}(:,:,3,3)=2*xi^2*E;
			tt{key}(:,:,4,3)=2*xi*eta*K;
			tt{key}(:,:,5,3)=2*xi*eta*L;
			tt{key}(:,:,2,4)=4*eta*K;
			tt{key}(:,:,4,4)=2*xi*E;
			tt{key}(:,:,2,5)=-4*eta*L;
			tt{key}(:,:,5,5)=2*xi*E;
		end
	end
end

tt=tt_matrix(tt); % @Bydlocode
return
end
