function mu2=calMu2(Pmat)

    [L,~,S,~] = size(Pmat);
    mu2 = size(S^2,L);
    for l = 1:L
        for z1 = 1:S
            for z2 = 1:S
                ss = (z1-1)*S + z2;
                if l~=1
                    v = L - sum(Pmat(l-1,1:(l-1),z1,1)) - sum(Pmat(l,l:end,z2,1));
                else
                    v = L - sum(Pmat(l,l:end,z2,1));
                end
                mu2(ss,l) = v;
            end
        end
    end
end