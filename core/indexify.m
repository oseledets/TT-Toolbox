function idx_x=indexify(xval,list_of_coresizes)
%this function takes an integer multiindex (in range 1..N-1) and unbundles it 
%into a set of smaller indices (in ranges: 1..coresize) according to core sizes
%If all core sizes are equal to 2, this function is
%equivalent to de2bi
xval=xval-1;
    for n=1:numel(xval)
        for i=1:size(list_of_coresizes,2)
            size_x=list_of_coresizes(n,i);
            idx_x(n,i)=mod(xval(n),size_x);
            xval(n)=xval(n)-idx_x(n,i);
            xval(n)=xval(n)/size_x;
        end
    end
idx_x=idx_x+1;
end