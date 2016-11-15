function y=purify(tt)
   %purify v2 unstable version
   %takes a ttv and checks if some cores are 1x1 in their i-index
   %if there are some,it contrapts them with adjacent cores, and removes those one-number index cores 

   G=core2cell(tt);
      
   i=1;
   while i<=numel(G) %it looks like numel(G) is not calculated on each step
        if size(G{i},2)==1
            if i<numel(G)
                %merge into right
                G1=G{i};
                G2=G{i+1};
                s1=size(G1,1);
                s2=size(G2,2);
                s3=size(G2,3);
                G1=reshape(G1,[s1*size(G1,2),size(G1,3)]);
                G2=reshape(G2,[size(G2,1),s2*s3]);
                GG=G1*G2;
                GG=reshape(GG,[s1,s2,s3]);
                G{i+1}=GG;
                G(i)=[];
            else %merge into left
                if numel(G)==1
                    y=tt;
                else
                    G1=G{i-1};
                    G2=G{i};
                    s1=size(G1,1);
                    s2=size(G2,2);
                    s3=size(G2,3);
                    keyboard
                    G1=reshape(G1,[s1*size(G1,2),size(G1,3)]);
                    G2=reshape(G2,[size(G2,1),s2*s3]);
                    GG=G1*G2;
                    GG=reshape(GG,[s1,s2,s3]);
                    G{i-1}=GG;
                    G(i)=[];
                end
            end
        else
                i=i+1;
        end
   end
        y=cell2core(tt_tensor,G);
end