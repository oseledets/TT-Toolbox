function [y]=tt_mat_full_vec(tt,x)
%Multiplication of the TT-matrix by full vector in TT1.0 format
%   [Y]=TT_MAT_FULL_VEC(TT,X) Multiplies full TT-matrix TT  by full 
%   vector X. Please avoid the usage: use "*" operator from the
%   object-oriented TT2 format. Will be removed in future releases.
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
d=size(tt,1);
n=size(x,1);
sz=zeros(1,d);
for i=1:d
    sz(i)=size(tt{i},1);
end
x=reshape(x,sz); 
%Convolve x over  mode-1 with tt{1}
n1=sz(1);
n2=n/n1;
core=tt{1};
ns1=size(core,1);
ns2=size(core,2);
r1=size(core,3);
core=permute(core,[1,3,2]); %core is i_1 \alpha_1 j_1
core=reshape(core,[ns1*r1,ns2]);
x=(core)*reshape(x,[n1,n2]); %x is i_1 \alpha_1,j_2,...,j_d
x=reshape(x, [n1, n2*r1]);
x=x.'; % x is \alpha_1,j_2, ...., j_d, i_1
msize=(n2*n1)/sz(2);
%msize=n2/sz(2);
%keyboard;
for k=2:d-1
    % x is alpha_{k-1},j_k,...,j_d,i_1,...,i_{k-1} convolve over alpha_{j-1} j_k
    core=tt{k}; %core is i_k j_k \alpha_{k-1} \alpha_k
    is=size(core,1); js=size(core,2); r1=size(core,3); r2=size(core,4);
    core=permute(core,[3,2,1,4]); %core is \alpha_{k-1} j_k i_k \alpha_k 
    core=reshape(core,[js*r1,is*r2]);
    x=reshape(x,[js*r1,msize]); 
    %keyboard;
    x=core.'*x; %x is i_k \alpha_k j_{k+1} ... j_d i_1,i_2,...,i_{k-1} 
    x=reshape(x,[is,msize*r2]);
    x=x.';
    msize=(is*msize)/js;
    %keyboard;
    %msize=
end
%x is \alpha_{d-1} j_d i_1 .... i_{d-1}
core=tt{d}; ns1=size(core,1); ns2=size(core,2); r=size(core,3);
core=permute(core,[3,2,1]); %core is \alpha_{d-1} j_d i_d  
core=reshape(core,[ns2*r,ns1]); 
x=reshape(x,[r*ns2,msize]);
x=core.'*x; %x is i_d i_1 ... i_{d-1}
x=reshape(x,[ns1,msize]);
x=x.';
y=reshape(x,[n,1]);

return
end
