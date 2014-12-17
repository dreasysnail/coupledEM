function [pCondStates] = hmmCondDecode(numStates,seq,pStates,pStates2)
%L of T
%Condition on 2 consequential states. P(z_i|z_l,x,theta)
%modified 2014-03-31 to speed up the calculation
L = length(seq);

pCondStates = zeros(L,L,numStates,numStates);

%Initiation
for i=1:L
    pCondStates(i,i,:,:) = eye(numStates,numStates);
end




%for speed consideration
%put the idx in last
repm = zeros(numStates,numStates,L);
for i = 1:L
    repm(:,:,i) = repmat(pStates(:,i),1,numStates);
end





%forward propagation
for l=1:L
    temp = eye(numStates,numStates);
    for i=l+1:L
        %T = (reshape(pStates2(:,i-1),numStates,numStates)./repmat(pStates(:,i-1)',2,1))';
        T = reshape(pStates2(:,i),numStates,numStates)'./repm(:,:,i-1);
        %divided by zero in pstates
        T(isnan(T))=0;
        T(isinf(T))=1;
        temp = temp*T;
        pCondStates(l,i,:,:) = temp;
        %assert(~any(any(isnan(temp))));
    end
end



%Backward propagation
for l=1:L
    temp = eye(numStates,numStates);
    for i=l-1:-1:1
        %T = (reshape(pStates2(:,i),numStates,numStates)./repmat(pStates(:,i+1),1,2));
        T = reshape(pStates2(:,i+1),numStates,numStates)./repm(:,:,i+1);
        %divided by zero in pstates
        T(isnan(T))=0;
        T(isinf(T))=1;
        temp = temp*T;
        pCondStates(l,i,:,:) = temp;
        %assert(~any(any(isnan(temp))));
    end
end

end








