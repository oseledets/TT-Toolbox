function [tt_mat]=nd_lapl2(lapl,d,d0,eps1)
%[TT_MAT]=ND_LAPL2(LAPL,D,D0,EPS1)
%D0-dimensional Laplace operator with N=2^D grid
%LAPL --- in TT format
%of size N=2^D, D is given
%D0 --- the dimension of the Laplace operator
%EPS1 --- truncation parameter for 1D Laplace (for the simplest 
%case the result is independent of it)
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

e=tt_eye(size(lapl{1},1),d);
e=tt_mat_to_vec(e);
lo=lapl;
lo=tt_mat_to_vec(lo);
tt_mat=[];
for k=1:d0
   tt_arr=cell(d0,1);
   for s=1:d0
      if ( s == k ) 
        tt_arr{s}=lo;
      else
        tt_arr{s}=e;
      end
   end
   %keyboard; 
   w=tt_mkron(tt_arr);
  if ( isempty(tt_mat) )
    tt_mat=w;
  else
    tt_mat=tt_add(tt_mat,w);
  end
  tt_mat=tt_compr2(tt_mat,eps1);
end
return
end

