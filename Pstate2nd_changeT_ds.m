function pState2 = Pstate2nd_changeT_ds(numStates,seq,fs,bs,scale,T,E,changeT)
    % L of T
    %zl,zl+1|x,theta
    seqLength = length(seq);
    seq = [0 seq];
    pState2 = zeros(numStates^2,seqLength);
    logf = log(fs);
    logb = log(bs);
    logE = log(E);
    len = (numStates-1)/2;
    TR = wrap_T(T(:,1),len);
    %complex number come from: TR<0
    for i = 1:seqLength
        if changeT
            TR(1,1) = 1-sum(T(:,i));
            TR(1,2) = T(1,i);
            TR(1,2+len) = T(2,i);
        end
        for k = 1:numStates
            for l = 1:numStates
                if TR(k,l)~=0              
                    pState2((k-1)*numStates+l,i) = exp( logf(k,i) + log(TR(k,l)) + logE(l,seq(i+1)) + logb(l,i+1))./scale(i+1);
                end
            end
        end
    end
end
