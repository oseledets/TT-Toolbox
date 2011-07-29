function [sub]=ind2sub2(siz,ind)
% function [sub]=ind2sub2(siz,ind)
% WORKING, bljat', ind2sub
% Splits integer index IND to a multiindex SUB according to dimensions SIZ

d = max(size(siz));
sub = zeros(1,d);
curnum = ind-1;
for i=d:-1:2
    sub(i) = floor(curnum/prod(siz(1:i-1)))+1;
    curnum = mod(curnum, prod(siz(1:i-1)));
end;
sub(1)=curnum+1;

end