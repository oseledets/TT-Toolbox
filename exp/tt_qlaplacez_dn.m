function [tt]=tt_qlaplacez_dn(d)

% returns a rank-8 QTT decomposition of
% Delta_{1}^{-1} \otimes \Id_{2} \ otimes \ldots \otimes \Id_{D} + \ldots +
%  + \Id_{1} \ otimes \ldots \otimes \Id_{D-1} \otimes Delta_{D}^{-1},
% Delta_{k} being a discretization of Laplace operator on 2^{d(k)} points
% uniform grid,
% Dirichlet boundary conditions being imposed
%
% D=size(d,2) must be >= 1
%
% September 3, 2010
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
J=zeros(2);
J(1,2)=1;
I2=zeros(2);
I2(2,2)=1;
E=ones(2);

if (D == 1)
	for key=1 : d
		tt{key}=eye(2);
	end
else

	key=0;
	for k=1 : D
		for kappa=1 : d(k)
			key=key+1;
			if (kappa == 1)
				if (k == 1)
					tt{key}=zeros(2,2,5);
					tt{key}(:,:,1)=I+I2+J+J';
					tt{key}(:,:,2)=2*E;
					tt{key}(:,:,3)=I2+J'+E;
					tt{key}(:,:,4)=I2+J+E;
					tt{key}(:,:,5)=I;	
				elseif (k == D)
					tt{key}=zeros(2,2,2,4);
					tt{key}(:,:,1,1)=I;
					tt{key}(:,:,2,1)=I+I2+J+J';
					tt{key}(:,:,2,2)=2*E;
					tt{key}(:,:,2,3)=I2+J'+E;
					tt{key}(:,:,2,4)=I2+J+E;
				else
					tt{key}=zeros(2,2,2,8);
					tt{key}(:,:,1,1)=I+I2+J+J';
					tt{key}(:,:,1,2)=2*E;
					tt{key}(:,:,1,3)=I2+J'+E;
					tt{key}(:,:,1,4)=I2+J+E;
					tt{key}(:,:,1,5)=I;
					tt{key}(:,:,2,5)=I+I2+J+J';
					tt{key}(:,:,2,6)=2*E;
					tt{key}(:,:,2,7)=I2+J'+E;
					tt{key}(:,:,2,8)=I2+J+E;
				end
			elseif (kappa == d(k))
				if (k == D)
					tt{key}=zeros(2,2,4);
					tt{key}(:,:,1)=I;
					tt{key}(:,:,2)=I2;
					tt{key}(:,:,3)=J;
					tt{key}(:,:,4)=J';
				elseif (k == 1)
					tt{key}=zeros(2,2,5,2);
					tt{key}(:,:,1,1)=I;
					tt{key}(:,:,2,1)=I2;
					tt{key}(:,:,3,1)=J;
					tt{key}(:,:,4,1)=J';
					tt{key}(:,:,5,2)=I;
				else
					tt{key}=zeros(2,2,8,2);
					tt{key}(:,:,1,1)=I;
					tt{key}(:,:,2,1)=I2;
					tt{key}(:,:,3,1)=J;
					tt{key}(:,:,4,1)=J';
					tt{key}(:,:,5,2)=I;
					tt{key}(:,:,6,2)=I2;
					tt{key}(:,:,7,2)=J;
					tt{key}(:,:,8,2)=J';
				end
			else
				if (k == D)
					tt{key}=zeros(2,2,4,4);
					tt{key}(:,:,1,1)=I;
					tt{key}(:,:,2,2)=2*E;
					tt{key}(:,:,3,3)=E;
					tt{key}(:,:,4,4)=E;
					tt{key}(:,:,2,1)=I2;
					tt{key}(:,:,3,1)=J;
					tt{key}(:,:,4,1)=J';
					tt{key}(:,:,2,3)=I2+J';
					tt{key}(:,:,2,4)=I2+J;
				elseif (k == 1)
					tt{key}=zeros(2,2,5,5);
					tt{key}(:,:,1,1)=I;
					tt{key}(:,:,2,2)=2*E;
					tt{key}(:,:,3,3)=E;
					tt{key}(:,:,4,4)=E;
					tt{key}(:,:,2,1)=I2;
					tt{key}(:,:,3,1)=J;
					tt{key}(:,:,4,1)=J';				
					tt{key}(:,:,2,3)=I2+J';
					tt{key}(:,:,2,4)=I2+J;
					tt{key}(:,:,5,5)=I;
				else
					tt{key}=zeros(2,2,8,8);
					tt{key}(:,:,1,1)=I;
					tt{key}(:,:,2,2)=2*E;
					tt{key}(:,:,3,3)=E;
					tt{key}(:,:,4,4)=E;
					tt{key}(:,:,2,1)=I2;
					tt{key}(:,:,3,1)=J;
					tt{key}(:,:,4,1)=J';
					tt{key}(:,:,2,3)=I2+J';
					tt{key}(:,:,2,4)=I2+J;
					tt{key}(:,:,5,5)=I;
					tt{key}(:,:,6,6)=2*E;
					tt{key}(:,:,7,7)=E;
					tt{key}(:,:,8,8)=E;
					tt{key}(:,:,6,5)=I2;
					tt{key}(:,:,7,5)=J;
					tt{key}(:,:,8,5)=J';
					tt{key}(:,:,6,7)=I2+J';
					tt{key}(:,:,6,8)=I2+J;
				end
			end
		end
	end
end

tt=tt_matrix(tt); % @Bydlocode
return
end
