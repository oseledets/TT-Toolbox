function newten= tenten_conv(ten1,d1, ind1, ten2, d2, ind2)
%  NEWTEN = TENTEN_CONV(TEN1,D1, IND1, TEN2, D2, IND2)
%  Convolution of two tensors
%  TEN1 with SIZE1 dimensions
%  TEN2 with SIZE2 dimensions
%  The sequence of tensors TEN1, TEN2 determines the sequence of indeces in
%  the resulting tensor
% 
%  The sets IND1 and IND2 include the convolution indeces
%  Output: the resulting tensor 
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


 %sten1=size(ten1);
 %sten2=size(ten2);

size1=d1;
size2=d2;
%keyboard;
col_conv_ind1=max(size(ind1));

ten_ind1=1:size1;
ten_ind1(ind1)=[];
ten_ind1=[ten_ind1, ind1];

% permute indices in ten1
ten1=permute(ten1, ten_ind1);

%mod_size1=zeros(1,size1-col_conv_ind1);
%full_size1=1;
sz1=size(ten1);
mod_size1=sz1(1:size1-col_conv_ind1);
full_size1=prod(mod_size1);
%full_size1=prod(sz1(1:size1-col_conv_ind1));
%for i=1:size1-col_conv_ind1
%    mod_size1(i)=size(ten1,i);
%    full_size1=full_size1*mod_size1(i);
%end

%conv_size=1;
conv_size=numel(ten1)/full_size1;
%for i=1:col_conv_ind1
%    conv_size=conv_size*size(ten1,size1+1-i);
%end

ten1=reshape(ten1, full_size1, conv_size);

% Now ten1 is a matrix   full_size1 x conv_size

col_conv_ind2=max(size(ind2));

ten_ind2=1:size2;
ten_ind2(ind2)=[];
ten_ind2=[ind2, ten_ind2];

% permute indices in ten2
ten2=permute(ten2, ten_ind2);
%size(ten2)

%mod_size2=zeros(1,size2-col_conv_ind2);
%full_size2=1;
sz2=size(ten2);
mod_size2=sz2(col_conv_ind2+1:size2);
full_size2=prod(mod_size2);
%for i=1:size2-col_conv_ind2
%    mod_size2(i)=size(ten2,col_conv_ind2+i);
%    full_size2=full_size2*mod_size2(i);
%end


ten2=reshape(ten2, conv_size, full_size2);

% Now ten2 is a matrix   conv_size  x  full_size2

newten=ten1*ten2;
newten=reshape(newten, [mod_size1, mod_size2]);

%snew=size(newten);

return
end
