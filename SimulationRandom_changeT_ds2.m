%simulating
function [seqs,Intensity,states,B,alpha,beta,sigma,bound,T,mu] = SimulationRandom_changeT_ds2(n,len,changeT)
%using E# as linker
    l        = 20;
    [alignmentSetMatrix,~] = readAlignmentSetMatrixFromFile('PWM/inputFiles/input.txt');
    [PFM, ~] = PWM( alignmentSetMatrix );
    index    = [3,1,4,2];
    PFM      = PFM(index,1:len)';
    %substitute
    %PFM = [0.271530000000000,0.111186000000000,0.145423000000000,0.471861000000000;0.289388000000000,0.117744000000000,0.137428000000000,0.455440000000000;0.265435000000000,0.114924000000000,0.131731000000000,0.487910000000000;0.125645000000000,0.555972000000000,0.0846770000000000,0.233706000000000;0.303948000000000,0.159204000000000,0.331508000000000,0.205341000000000;0.302736000000000,0.261673000000000,0.225970000000000,0.209621000000000];
    PFM_r    = PFM(len:-1:1,4:-1:1);
    BG       = [.25,.25,.25,.25];
    B        = [BG;PFM;PFM_r];
    for i = 1:n
        seqs(i,:) = randsample(1:4,l,true,BG);
    end
%%use true B
    [T,~,~] = hmmtrain_fixE(seqs,wrap_T([0.0011,0.0011],6),B,'Verbose',true,'Maxiterations',5);
     
    if changeT 
        T = betapdf(0.001:1/l:0.999-1/l,6,2);
        T = T/max(T)/4;
    else
        T = 1-T(1,1);
    end
    bound = [];
    %s1 = [];
    %s2 = [];
    T = [T/2;T/2];
    
    alpha = 10;
    beta = 4;
    sigma = 1;

    for i = 1:n
        [pStates,logPseq, ~,~, ~] = hmmdecode_changeT_ds(seqs(i,:),T,B,changeT);
        bound = [bound;1-pStates(1,:)];
        %s1 = [s1,sum(pStates(2:len+1,:),1)];
        %s2 = [s2,sum(pStates(len+2:end,:),1)];
        [currentState, logP] = hmmviterbi_topK(seqs(i,:),wrap_T(T,len),B,30);
        pp = exp(logP-logPseq);
        if sum(pp)<0.95
            fprintf('Warning! Viterbi paths contain only %d prob of %dth seq\n',sum(pp),i);
        end
        pp = pp./sum(pp);
        boundNum = sum(currentState~=1,1);
        mu(i) = randsample(boundNum,1,true,pp);
        Intensity(i) = alpha + mu(i)*beta + normrnd(0,sigma) ;
        % p(#|x,y)\propto p(#|x)p(y|#)
    end
    
    states = bound; 

    Intensity = Intensity';
    %??
%     if ~changeT
%         [guessTR,~,~] = hmmtrain_fixE(seqs,wrap_T(T,len),B);
%     else
%         [guessTR,~,~] = hmmtrain_changeT_fixE(seqs,T,B);
%     end
   
    
    
   