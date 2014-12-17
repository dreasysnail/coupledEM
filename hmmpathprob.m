function p=hmmpathprob(path,seq,TR,E,logpseq)
%path \in {1,2,...}
    l = length(seq);
    p = log(TR(1,path(1)));
    for i = 1:l-1
        p = p+log(TR(path(i),path(i+1)));
        p = p+log(E(path(i),seq(i)));
    end
    p = p+log(E(path(l),seq(l)));
    p = p - logpseq;
    p = exp(p);