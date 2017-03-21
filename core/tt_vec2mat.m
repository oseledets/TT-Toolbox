function ttm=tt_vec2mat(ttv)
   
    %    ttv is a QTT tt_vector(tt_tensor) of a full lenght N^2
    %    ttm is a QTT tt_matrix equivalent to simple reshape of this vector
    %    into matrix

    % TT-Toolbox 2.2.1, 2009-2016
    %
    %This is TT Toolbox, written by Ivan Oseledets et al.
    %Institute of Numerical Mathematics, Moscow, Russia
    %webpage: http://spring.inm.ras.ru/osel
    %
    %For all questions, bugs and suggestions please mail
    %ivan.oseledets@gmail.com
    %---------------------------
    %Last tweaks by Alexey Boyko of Oseledets' group
    %alexey.boyko@skolkovotech.ru    

%     d=ttv.d/2; 
    modes  = ttv.n(1:2:ttv.d);
    modes2 = ttv.n(2:2:ttv.d);

%     list_for_standard_permute=[1:d;d+1:d+d];
%     list_for_standard_permute=list_for_standard_permute(:)';
    ttv_readyformat=reshape(ttv,modes.*modes2);
    ttm=tt_matrix(ttv_readyformat,modes,modes2);
end