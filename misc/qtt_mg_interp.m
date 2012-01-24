function [m]=qtt_mg_interp(d)
%  function [m]=qtt_mg_interp(d)
% Generates the dumb multigrid interpolator (tt_matrix)
% 1,	0
% 0.5,	0.5
% 0,	1
% 0,	0.5

prol0 = [1 0; 0.5 0.5; 0 1; 0 0.5];
h0 = [0 0; 0 0; 0 0; 0.5 0];
I = eye(2);
I12 = [0 1; 0 0];
I21 = [0 0; 1 0];

m = cell(d,1);

m{1}=zeros(4,2,2);
m{1}(:,:,1)=prol0; m{1}(:,:,2)=h0;
for i=2:d-1
    m{i}=zeros(2,2,2,2);
    m{i}(:,:,1,1)=I;
    m{i}(:,:,2,1)=I12; m{i}(:,:,2,2)=I21;
end; 
m{d} = zeros(2,2,2);
m{d}(:,:,1)=I; m{d}(:,:,2)=I12;

m = tt_matrix(m);

end