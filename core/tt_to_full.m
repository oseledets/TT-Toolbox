function [f]=tt_to_full(tt)
%[F]=TT_TO_FULL(TT)
%Compute full array F from its tensor-train representation TT
%Inverse of full_to_tt
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
%d=size(tt,1);



d=size(tt,1);
f=tt{d};
bs=size(f,1);
f=f.';
for i=d-1:-1:2
    core=tt{i}; 
    ncur=size(core,1);
    r2=size(core,2);
    r3=size(core,3);
    %core is [ncur,r2,r3], f is r3xbs
    core=reshape(core,[ncur*r2,r3]);
    f=core*f; %f is ncur*r2*bs, should be r2*ncur*bs
    f=reshape(f,[ncur,r2,bs]);
    f=permute(f,[2,1,3]);
    f=reshape(f,[r2,ncur*bs]);
    bs=ncur*bs;    
end
f=tt{1}*f;
sz=zeros(1,d);
for i=1:d
    sz(i)=size(tt{i},1);
end
f=reshape(f,sz);
return
end
