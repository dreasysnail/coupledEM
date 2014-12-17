function [kmer,avg] = medianFreq (seqs,intensity,len)
    [numSeqs,seqLength] = size(seqs);
    dic = containers.Map;
    count = containers.Map;
    for i = 1:numSeqs
        for l = 1:seqLength-len
            mer = mat2str(seqs(i,l:(l+len-1)));
            if isKey(dic,mer)             
                dic(mer) = dic(mer) + intensity(i);
                count(mer) = count(mer)+1;
            else
                dic(mer) = intensity(i);
                count(mer) = 1;
            end
        end
    end
    kmertmp = keys(dic);
    kmer = zeros(length(kmertmp),len);
    avg = zeros(length(kmertmp),1);
    for i = 1:length(kmertmp)
        kmer(i,:) = str2num(str2mat(kmertmp(i)));
        avg(i) = dic(kmertmp{i})/count(kmertmp{i});
    end
    [avg,idx] = sort(avg,'descend');
    kmer = kmer(idx,:);
end
