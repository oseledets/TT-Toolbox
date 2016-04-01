function y1=tt_submatrix(ttm,list_of_idx,list_of_vals)
%% TT-Toolbox 2.2.2, 2009-2016
% list_of_idx should be in strict ascending order!
% list of vals is 2-by-d array containing values in range from 1 to coresize{i}

%This function freezes some (or all) cores to some particular values
% from appropriate values from 2-by-numel(list_of_idx) variable list_of_vals
% assymetric submatrix extraction (for example, extract a rectangular
% matrix from a square matrix) is not yet supported

%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel

%for questions on this particular code email Alexey Boyko,
% Skolkovo Institute of Science and Technology,
% a.boyko@skoltech.ru or alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2016

G=core2cell(ttm);

Nidx=numel(list_of_idx);

for k=Nidx:-1:1
    idx = list_of_idx(k);
    val = list_of_vals(:,k);
    G{idx}=G{idx}(:,val(1),val(2),:);

    if numel(G)>1
        BackConvolve = (idx>1); 
%     if numel(G{idx})~=1
        if BackConvolve
            %here is the scenario fo convolve with previous index   
            G{idx}=permute(G{idx},[1 4 2 3]);
            sz=size(G{idx-1});
            sz=[sz ones(1,4-numel(sz))];

            G{idx-1}=reshape(G{idx-1},[sz(1)*sz(2)*sz(3),sz(4)]);
            G{idx-1}=G{idx-1}*G{idx};
            G=G(setdiff(1:numel(G), idx) ); %we keep all core numbers except for the frozen one (idx)
            G{idx-1}=reshape(G{idx-1},[sz(1),sz(2),sz(3),size(G{idx-1},2)]);
        else
            %here is the scenario fo convolve with forward index           
            G{idx}=permute(G{idx},[1 4 2 3]);
            sz=size(G{idx+1});
            sz=[sz ones(1,4-numel(sz))];
            G{idx+1}=reshape(G{idx+1},[sz(1),sz(2)*sz(3)*sz(4)]);
            G{idx+1}=G{idx}*G{idx+1};
            G=G(setdiff(1:numel(G), idx) );
            %we keep all core numbers except for the frozen one (idx)
            G{idx}=reshape(G{idx},[size(G{idx},1),sz(2),sz(3),sz(4)]);
        end
    end   
end

if numel(G) ==1 && numel(G{:})==1
    y1=G{:};
else
    y1=cell2core(tt_matrix,G);
end

end


