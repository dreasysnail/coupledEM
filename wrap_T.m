function A = wrap_T(T,len)
%T =  ( T(1),T(2) )
    A = zeros(2*len+1,2*len+1);
    %Modification of A
    A(1,1) = 1-sum(T);
    A(1,2) = T(1);
    A(1,len+2) = T(2);
    for j = 2:len
        A(j,j+1) = 1;
    end
    for j = len+2:2*len
        A(j,j+1) = 1;
    end
    A(len+1,1) = 1;
    A(2*len+1,1) = 1;
end

% function A = wrap_T(T,len)
%     A        = reshape([[1-T,T],zeros(1,(len+1)^2-2)],len+1,len+1)';
%     %Modification of A
%     for j = 2:len
%         A(j,j+1) = 1;
%     end
%     A(len+1,1) = 1;
% end
