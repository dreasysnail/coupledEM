function [id,val] = topK_in_matrix(mat,k)
    m = size(mat,1);
    %[val,maxIndex] = maxk(mat(:),k); 
    [sortedValues,sortIndex] = sort(mat(:),'descend'); 
    val = sortedValues(1:k);
    maxIndex = sortIndex(1:k); 
    id(:,1) = rem(maxIndex-1,m) + 1;
    id(:,2) = (maxIndex-id(:,1))/m + 1;   
%     id2 = [mod(maxIndex,m),floor(maxIndex/m)+1];
%     idx = (id2(:,1)==0);
%     id2(idx,:) = [id2(idx,1) + m,id2(idx,2) - 1];
end
    