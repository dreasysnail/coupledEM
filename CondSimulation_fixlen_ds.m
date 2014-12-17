%simulating
function [seqs,Intensity,states,B,alpha,beta,sigma,bound,T] = CondSimulation_fixlen_ds(n,len,changeT)
%L of T    
%[seqs,states,Intensity,Beta] = simulation_fixlen(n,len)
    addpath(genpath(pwd));
    l        = 20;
    t        = 0.25;

    %Modification of A
    
    BG       = [.25,.25,.25,.25];
    [alignmentSetMatrix,~] = readAlignmentSetMatrixFromFile('PWM/inputFiles/input.txt');
    [PFM, ~] = PWM( alignmentSetMatrix );
    
    %% plot3 PWM

    %[PFMsorted, ~, Idx_row] = sortColumnsByBaseOccurenceFreq(PFM(:,1:len));

    %sizeAndPositionDNA_bases (PFMsorted, Idx_row);


    index    = [3,1,4,2];
    PFM      = PFM(index,1:len)';
    PFM_r    = PFM(len:-1:1,4:-1:1);
    B        = [BG;PFM;PFM_r];
    
    %B = [0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000;0.271530000000000,0.111186000000000,0.145423000000000,0.471861000000000;0.289388000000000,0.117744000000000,0.137428000000000,0.455440000000000;0.265435000000000,0.114924000000000,0.131731000000000,0.487910000000000;0.125645000000000,0.555972000000000,0.0846770000000000,0.233706000000000;0.303948000000000,0.159204000000000,0.331508000000000,0.205341000000000;0.302736000000000,0.261673000000000,0.225970000000000,0.209621000000000;0.209621000000000,0.225970000000000,0.261673000000000,0.302736000000000;0.205341000000000,0.331508000000000,0.159204000000000,0.303948000000000;0.233706000000000,0.0846770000000000,0.555972000000000,0.125645000000000;0.487910000000000,0.131731000000000,0.114924000000000,0.265435000000000;0.455440000000000,0.137428000000000,0.117744000000000,0.289388000000000;0.471861000000000,0.145423000000000,0.111186000000000,0.271530000000000];


    
    
   
    seqs     = zeros(n,l);
    states   = zeros(n,l);
    
    alpha = 10;
    beta = 4;
    sigma = 5;
    %last T is useless   length(T)=l
    if changeT
        T        = betapdf(0.001:1/l:0.999,6,2);
        T = T/max(T)*t+0.1;
    else
        T        = ones(l,1)*t+0.1;
    end

    for i = 1:n
        [seqs(i,:),states(i,:)] = hmmgenerate_changeT_ds(l,wrap_T([T(1)/2,T(1)/2],len),B,T) ;
        states(i,:) = states(i,:)-1 ;
        bound(i,:) = states(i,:)>0 +0;
        % or the ensemble bound number?
        Intensity(i) = alpha + sum(bound(i,:))*beta + normrnd(0,sigma) ;
    end
    if ~changeT
        T = t+0.1;
    end
    Intensity = Intensity';










   
