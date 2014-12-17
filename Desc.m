function rs=Desc(set,k)
%Desc product
    if k==1
        rs=cell2mat(set');
        return;
    end
    rs=[];
    for i=1:length(set)
        mat=Desc(set,k-1);
        mat=[set{i}.*ones(size(mat,1),1),mat];
        rs=[rs;mat];
    end
end