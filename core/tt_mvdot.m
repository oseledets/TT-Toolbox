function [res] = tt_mvdot(a_tt,x_tt,y_tt)
%Fast computation of the bilinear form (Ax,y) in the TT-format
%   [RES]=TT_MVDOT(A_TT,X_TT,Y_TT) Computes scalar product(AX,Y), where 
%   A_TT is a matrix, X_TT is a vector, Y_TT  is a vector
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
d=size(a_tt,1);
%The c
core=a_tt{1}; matx=x_tt{1}; maty=y_tt{1};
n1=size(core,1); n2=size(core,2);
rc=size(core,3);
rx=size(matx,2); ry=size(maty,2);
core=maty'*reshape(core,[n1,n2*rc]); %core is ryxn2xrc
core=reshape(permute(reshape(core,[ry,n2,rc]),[3,1,2]),[rc*ry,n2]);
phi=permute(reshape(core*matx,[rc,ry,rx]),[1,3,2]); %core is rcxryxrx-> rcxrx*ry

for i=2:d-1
  core=a_tt{i}; core_x = x_tt{i}; core_y=y_tt{i};
  %Decipher sizes
 n1=size(core,1); n2=size(core,2); rc1=size(core,3); rc2=size(core,4);
 rx1=size(core_x,2); rx2=size(core_x,3); ry1=size(core_y,2); ry2=size(core_y,3);
 %Convolve phi & core over the first index
 %core is n1xn2xrc1xrc2 phi is rc1x rx1 x ry1
  core=reshape(permute(core,[3,1,2,4,5]),[rc1,n1*n2*rc2]);
  phi=reshape(phi,[rc1,rx1*ry1])'*core; %Is rx1 x ry1 x n1 x n2 x rc2
  %Convolve over n2 & rx1 
  phi=permute(reshape(phi,[rx1,ry1,n1,n2,rc2]),[4,1,2,3,5]); 
  %n2 x rx1 x ry1 x n1 x rc2;
  phi=reshape(core_x,[n2*rx1,rx2])'*reshape(phi,[n2*rx1,ry1*n1*rc2]);
  %phi is now rx2 x ry1 x n1 x rc2
  %Convolve over n1 & ry1 
  phi=reshape(phi,[rx2,ry1*n1,rc2]); phi=permute(phi,[3,1,2]);
  %keyboard;
  phi=reshape(phi,[rc2*rx2,ry1*n1])*reshape(permute(core_y,[2,1,3]),[ry1*n1,ry2]);
  phi=reshape(phi,[rc2,rx2,ry2]);
end
  core=a_tt{d}; mat_x=x_tt{d}; mat_y=y_tt{d};
  n1=size(core,1); n2=size(core,2); rc=size(core,3);
  rx=size(mat_x,2); ry=size(mat_y,2);
  %Convolve phi & core over aux index
  %core is n1 x n2 x rc
  phi=reshape(core,[n1*n2,rc])*reshape(phi,[rc,rx*ry]);
  %phi is n1 x n2 x rx x ry
  %Convolve over n2 & rx
  phi=permute(reshape(phi,[n1,n2*rx,ry]),[1,3,2]);
  phi=reshape(phi,[n1*ry,n2*rx])*reshape(mat_x,[n2*rx,1]);
  %phi is a column of length 
  res= phi'*reshape(mat_y,[n1*ry,1]);

return

end
