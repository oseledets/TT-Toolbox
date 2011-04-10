% test code for tt_qlaplacex_dd()
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
tt=tt_qlaplacex_dd(d);
tt=tt_mat_to_vec(tt);


full=nd_to_full(tt);

disp('inv computation... ');
Li=cell(D,1);
for k=1 : D
	L=2*eye(2^d(k));
	for i=1 : 2^d(k)-1
		L(i,i+1)=-1;
		L(i+1,i)=-1;
	end
	Li{k}=inv(L);
end
disp('OK');

disp('kron computation... ');
Z=Li{1};
for k=2 : D
	Z=kron(Z,Li{k});
end
disp('OK');

err=norm(full-Z,'fro');
fprintf('fro err = %e\n', err);
