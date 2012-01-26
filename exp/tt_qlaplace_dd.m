function [tt]=tt_qlaplace_dd(d)
% returns a rank-3,4...4 QTT decomposition of
% Delta_{1} \otimes \Id_{2} \ otimes \ldots \otimes \Id_{D} + \ldots +
%  + \Id_{1} \ otimes \ldots \otimes \Id_{D-1} \otimes Delta_{D},
% Delta_{k} being a discretization of Laplace operator on 2^d(k) points
% uniform grid,
% Dirichlet-Dirichlet boundary conditions being imposed
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

if (D == 1)
	for key=1 : d
		if (key == 1)
			tt{key}=zeros(2,2,3);
			tt{key}(:,:,1)=2*I-J-J';
			tt{key}(:,:,2)=-J;
			tt{key}(:,:,3)=-J';
		elseif (key == d)
			tt{key}=zeros(2,2,3);
			tt{key}(:,:,1)=I;
			tt{key}(:,:,2)=J';
			tt{key}(:,:,3)=J;
		else
			tt{key}=zeros(2,2,3,3);
			tt{key}(:,:,1,1)=I;
			tt{key}(:,:,2,2)=J;
			tt{key}(:,:,3,3)=J';
			tt{key}(:,:,2,1)=J';
			tt{key}(:,:,3,1)=J;
		end
	end
else

	key=0;
	for k=1 : D
		for kappa=1 : d(k)
			key=key+1;
			if (kappa == 1)
				if (k == 1)
					tt{key}=zeros(2,2,4);
					tt{key}(:,:,1)=2*I-J-J';
					tt{key}(:,:,2)=-J;
					tt{key}(:,:,3)=-J';
					tt{key}(:,:,4)=I;
				elseif (k == D)
					tt{key}=zeros(2,2,2,3);
					tt{key}(:,:,1,1)=2*I-J-J';
					tt{key}(:,:,1,2)=-J;
					tt{key}(:,:,1,3)=-J';
					tt{key}(:,:,2,1)=I;
				else
					tt{key}=zeros(2,2,2,4);
					tt{key}(:,:,1,1)=2*I-J-J';
					tt{key}(:,:,1,2)=-J;
					tt{key}(:,:,1,3)=-J';
					tt{key}(:,:,1,4)=I;
					tt{key}(:,:,2,1)=I;
				end
			elseif (kappa == d(k))
				if (k == D)
					tt{key}=zeros(2,2,3);
					tt{key}(:,:,1)=I;
					tt{key}(:,:,2)=J';
					tt{key}(:,:,3)=J;
				else
					tt{key}=zeros(2,2,4,2);
					tt{key}(:,:,4,1)=I;
					tt{key}(:,:,1,2)=I;
					tt{key}(:,:,2,2)=J';
					tt{key}(:,:,3,2)=J;
				end
			else
				if (k == D)
					tt{key}=zeros(2,2,3,3);
					tt{key}(:,:,1,1)=I;
					tt{key}(:,:,2,2)=J;
					tt{key}(:,:,3,3)=J';
					tt{key}(:,:,2,1)=J';
					tt{key}(:,:,3,1)=J;
				else
					tt{key}=zeros(2,2,4,4);
					tt{key}(:,:,1,1)=I;
					tt{key}(:,:,2,2)=J;
					tt{key}(:,:,3,3)=J';
					tt{key}(:,:,2,1)=J';
					tt{key}(:,:,3,1)=J;
					tt{key}(:,:,4,4)=I;
				end
			end
		end
	end
end
tt=tt_matrix(tt); % @Bydlocode
return
end