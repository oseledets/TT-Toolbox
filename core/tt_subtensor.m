function y1=tt_subtensor(ttv,list_of_idx,list_of_vals)
%% TT-Toolbox 2.2.1, 2009-2016
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel

%for questions on this particular code email Alexey Boyko, Skolkovo Institute of Science and Technology,
% a.boyko@skoltech.ru or alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2016
%%
% list_of_idx should be in strict ascending order!
% list of vals must contain 1 or 2 as values
G=core2cell(ttv);
% keyboard
Nidx=numel(list_of_idx);

for k=Nidx:-1:1
    idx = list_of_idx(k);
    val = list_of_vals(k);
    BackConvolve = (idx>1); 
    G{idx}=G{idx}(:,val,:);
    
%     if numel(G{idx})~=1
        if BackConvolve
            %here is the scenario fo convolve with previous index   
            G{idx}=permute(G{idx},[1 3 2]);
            sz=size(G{idx-1});
            if numel(sz)==2
                sz(3)=1;
            end
            G{idx-1}=reshape(G{idx-1},[sz(1)*sz(2),sz(3)]);
            G{idx-1}=G{idx-1}*G{idx};
            G=G(setdiff(1:numel(G), idx) ); %we keep all core numbers except for the frozen one (idx)
            G{idx-1}=reshape(G{idx-1},[sz(1),sz(2),size(G{idx-1},2)]);
        else
            %here is the scenario fo convolve with forward index           
            G{idx}=permute(G{idx},[1 3 2]);
            sz=size(G{idx+1});
            if numel(sz)==2
                sz(3)=1;
            end
            G{idx+1}=reshape(G{idx+1},[sz(1),sz(2)*sz(3)]);
            G{idx+1}=G{idx}*G{idx+1};
            G=G(setdiff(1:numel(G), idx) );
            %we keep all core numbers except for the frozen one (idx)
            G{idx}=reshape(G{idx},[size(G{idx},1),sz(2),sz(3)]);
        end
% end
        
end

if numel(G) ==1 && numel(G{:})==1
    y1=G{:};
else
    y1=cell2core(tt_tensor,G);
end


 

end


