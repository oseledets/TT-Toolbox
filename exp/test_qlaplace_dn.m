% test code for tt_qlaplace_dn()
%
% September 22, 2010
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

% d is the only parameter
d=[3,4,5];
%

D=size(d,2);
tt=tt_qlaplace_dn(d);
tt=tt_mat_to_vec(tt);


full=nd_to_full(tt);

Z=zeros(2^sum(d));
for k=1 : D
	L=2*eye(2^d(k));
	for i=1 : 2^d(k)-1
		L(i,i+1)=-1;
		L(i+1,i)=-1;
	end
	L(2^d(k),2^d(k))=1;
	for m=k-1 : -1: 1
		L=kron(eye(2^d(m)),L);
	end
	for m=k+1 : D
		L=kron(L,eye(2^d(m)));
	end
	Z=Z+L;
end
err=norm(full-Z,'fro');
fprintf('fro err = %e\n', err);
