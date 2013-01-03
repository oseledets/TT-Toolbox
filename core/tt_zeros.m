function [tt] = tt_zeros(n,varargin)
%Tensor of all zeros
%   [TT]=TT_ZEROS(N,D), computes the d-dimensional TT-tensor of all zeros 
%   with mode size equal to N
%
%   [TT]=TT_ZEROS(N), computes the TT-tensor of all zeros with mode size 
%   given in the vector N
%
%
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
if (numel(n) == 1)
  if (numel(varargin)>0)
      d=varargin{1}; 
  else
      d=1;
  end;
  n=n*ones(d,1);
else
  if(size(n,1)==1)
      n=n';
  end
  d=size(n,1);
end
tt=cell(d,1);
tt{1}=zeros([n(1,:),1]);
tt{d}=zeros([n(d,:),1]);
for i=2:d-1
   tt{i}=zeros([n(i,:),1,1]);
end
switch size(n,2)
    case 1
        tt=tt_tensor(tt); %Bydlocode
    case 2
        tt=tt_matrix(tt); %Bydlocode
    otherwise
        error('tt_zeros: N array has illegal size.\n');
end
return
end
