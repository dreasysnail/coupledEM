%final model
%T unchanged
len=6;
n=500;
changeT = false;
%[seqs,intensity,states.true,B.true,alpha.true,beta.true,sigma.true,bound,T.true,mu.true] = SimulationRandom_changeT_ds(n,len,changeT);
%[seqs,intensity,states.true,B.true,alpha.true,beta.true,sigma.true,bound,T.true,mu.true] = SimulationRandom_changeT_ds2(n,len,changeT);
%[guessTR,~,~] = hmmtrain_fixE(seqs,wrap_T(T.true,len),B.true,'Verbose',true);
[seqs,intensity,states.true,B.true,alpha.true,beta.true,sigma.true,bound,T.true]=CondSimulation_fixlen_ds(n,len,changeT);
%mu.true = sum(bound,2);
k = 1;
numStates = len*2+1;
numEmission =4;
T.guess = sum(T.true)*(1+k*rand(1,1));
% T.guess = T.true + unifrnd(0,1,size(T.true))*k;
% T.guess = T.guess./max(T.guess)*max(T.true);

%A.guess = A.guess./repmat(sum(A.guess,numStates),1,numStates);
B.guess = B.true + unifrnd(0,1,size(B.true))*k;
B.guess = B.guess./repmat(sum(B.guess,2),1,numEmission);
alpha.guess =alpha.true*(1+k*rand(1,1));
beta.guess = beta.true*(1+k*rand(1,1));
sigma.guess =sigma.true*(1+k*rand(1,1));
%A.dist = norm(A.guess - A.true,inf);
%B.dist = norm(B.guess - B.true,inf);

%[T.emA,B.emA,alpha.emA,beta.emA,sigma.emA,~] = hmmIntTrain_Motif(seqs,intensity,T.guess,B.guess,alpha.guess,beta.guess,sigma.guess,changeT);
[T.emA,B.emA,alpha.emA,beta.emA,sigma.emA,~] = hmmIntTrain_Motif_realDist(seqs,intensity,T.guess,B.guess,alpha.guess,beta.guess,sigma.guess,changeT);


%hmmIntTrain_Motif_realDist(seqs,intensity,sum(T.em),B.em,alpha.em,beta.em,sigma.em,changeT);
[T.smoothing,B.smoothing,~] = hmmtrain_changeT_ds(seqs,T.guess,B.guess,changeT);




%random plant motif

for count = 1:n
[pStates,logPseq, fs, bs, s] = hmmdecode(seqs(count,:),T.smoothing,B.smoothing);
expBound(count) = 20-sum(pStates(1,:)); 
end







%real data
%from bad init   from fixed T
%l=60;
k=0;
% T.guess = betapdf(0.001:1/l:0.999-1/l,6,2);
% T.guess = T.guess + unifrnd(0,1,size(T.guess))*k;
% T.guess = T.guess/max(T.guess)*0.25;
len = 6;
numStates=len*2+1;
numEmission =4;
changeT = false;
idx = randsample(1:10000,10000,false);
rearSeqs = seqs(idx,30:60);
rearlogInt = logInt(idx,:);
%T.truncatedguess = T.guess(40:59);
T.truncatedguess = 0.001;
B.guess = guessE + unifrnd(0,1,size(guessE))*k;
B.guess = B.guess./repmat(sum(B.guess,2),1,numEmission);

%%T.truncatedguess = 0.25;
[T.em,B.em,alpha.em,beta.em,sigma.em,logliks] = hmmIntTrain_Motif_realDist(rearSeqs,rearlogInt,T.truncatedguess,B.guess,7.3,0.1,0.3,changeT);
%[T.em,B.em,alpha.em,beta.em,sigma.em,logliks] = hmmIntTrain_Motif(rearSeqs,rearlogInt,oldGuessT,oldGuessE,alpha,beta,sigma,changeT);

% hmmIntTrain_Motif(seqs,intensity,T.em,B.em,alpha.em,beta.em,sigma.em,changeT);
[T.smoothing,B.smoothing,~] = hmmtrain_changeT_ds(rearSeqs,T.guess,B.guess);

% [kmer,avg] = medianFreq (rearSeqs,rearlogInt, len);

%%validate if it is linear
n =10000;
for i = 1:n
    [pStates,~, ~, ~, ~] = hmmdecode_changeT_ds(seqs(i,:),[T.truncatedguess/2;T.truncatedguess/2],B.guess,changeT);
    mu.true(i) = sum(1-pStates(1,:));
end
simpleRegression(logInt,mu.true)
scatter(mu.true,logInt)
scatter(mu.true,rearlogInt)
    
%rupa
BG = ones(1,4)*.25;
B.rupa = [BG;PFM;PFM(len:-1:1,4:-1:1)];


%scatter plot
mu.temp= [];
for i=1:n
    [ptemp,~, ~, ~, ~] = hmmdecode_changeT_ds(seqs(i,:),T.em,B.guess,changeT);
    mu.temp(i) = sum(1-ptemp(1,:));
end
    scatter(mu.temp,Intensity,'b.');
    hold on;
    xlabel('E(#|x,theta)');
    ylabel('log intensity');
    lsline;


%%heatmap
rl = {'b','m1','m2','m3','m4','m5','m6','c1','c2','c3','c4','c5','c6'};
cl= {'A','C','G','T'};
HeatMap(B.true,'RowLabels',rl,'ColumnLabels',cl,'Colormap','autumn','Symmetric',true);
HeatMap(B.em,'RowLabels',rl,'ColumnLabels',cl,'Colormap','autumn','Symmetric',true);


HeatMap(B.twostep,'RowLabels',rl,'ColumnLabels',cl,'Colormap','autumn','Symmetric',true);
HeatMap(B.rupa,'RowLabels',rl,'ColumnLabels',cl,'Colormap','autumn','Symmetric',true);

HeatMap(B.guess,'RowLabels',rl,'ColumnLabels',cl,'Colormap','autumn','Symmetric',true);




%posterior vs profile
        hold off;
        dbup
        hold on; bar(bound(4,:)+0,'g')
        dbdown
        %hold on; bar((seq==1)/10.0,'r')
        %hold on; bar((seq==3)/10.0,'b')
        plot(1-pStates(1,:),'b','LineWidth',2)
        hold on; plot(1-condLik(1,:)./(mean(condLik(2:end,:),1)+condLik(1,:)),'black--','LineWidth',2)
        %hold on; plot(condLik(1,:),'r--')
        hold on; plot(1-posterior(1,:),'r','LineWidth',2)
        %hold on; plot(1-posterior_tmp(1,:),'g','LineWidth',2)
        legend('True states','p(z|x)','Reweigth','p(z|x,y)')

        
%plot T
        hold off;
        plot(T.guess);
        hold on;
        plot(T.true,'b--');
        hold on;
        plot(T.smoothing,'g');
        plot(T.em,'r');
        xlabel('Position on Probe');
        ylabel('P(BG->Motif)');
        legend('Initial guess','Ground Truth','w/o Intensity','w/ Intensity',...
       'Location','NW');
   
%%Roc 
%%Plot2:Prob Profile
states.realbound = (states.real>0)+0;
for i=1:n
    [ptemp,~, ~, ~, ~] = hmmdecode_changeT_ds(seqs(i,:),T.true,B.true);
    states.trueParam(i,:) = 1-ptemp(1,:);
    [ptemp,~, ~, ~, ~] = hmmdecode_changeT_ds(seqs(i,:),T.em,B.em);
    states.em(i,:) = real(1-ptemp(1,:));
    [ptemp,~, ~, ~, ~] = hmmdecode_changeT_ds(seqs(i,:),T.smoothing,B.smoothing);
    states.smoothing(i,:) = real(1-ptemp(1,:));
end



[roc.tpr1,roc.fpr1,roc.thr] = ROC(states.realbound,states.trueParam);
[roc.tpr2,roc.fpr2,~] = ROC(states.realbound,states.em);
[roc.tpr3,roc.fpr3,~] = ROC(states.realbound,states.smoothing);
plot(roc.fpr1,roc.tpr1,'r',roc.fpr2,roc.tpr2,'b',roc.fpr3,roc.tpr3,'g');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
legend('True params','EM w/ Intensity','EM w/o Intensity',...
       'Location','NW');
   
   
%% Expected bound #


mu.em = sum(states.em,2); 
mu.smoothing = sum(states.smoothing,2); 
hold off;
scatter(intensity,mu.smoothing,'r');
hold on;
scatter(intensity,mu.em,'b');
xlabel('intensity');
ylabel('E(bound#)');



%%two step
 T.true = [T.true/2,T.true/2];
[T.twostep,B.twostep,~] = hmmtrain(seqs,wrap_T(T.true,6),B.guess,'Verbose',true,'Maxiterations',200);
T.twostep = real(T.twostep);
for i=1:n
    [ptemp,~, ~, ~, ~] = hmmdecode(seqs(i,:),T.twostep,B.twostep);
    states.twostep(i,:) = 1-ptemp(1,:);
end
mu.twostep = sum(states.twostep,2); 
[alpha.twostep,beta.twostep,sigma.twostep] = simpleRegression(intensity,mu.twostep);


   
   
%test distribution of p(#|x)
    l        = 20;
    %why A lead to qi tou

    dist =zeros(n,l+1);
    tempmu=zeros(n,1);
    tempvar=zeros(n,1);
    Paths = allPath(l-1,len);
    Paths = [zeros(size(Paths,1),1),Paths];
    T.true = ones(1,l)*T.guess;
    for i = 1:n
        seq = seqs(i,:);
        [pStates,logpseq, fs,bs, s] = hmmdecode_changeT_ds(seq,T.true,B.true,changeT);
%         tempmu(i) = l-sum(pStates(1,:)); 
%         pStates2 = Pstate2nd_changeT_ds(numStates,seq,fs, bs,s,T.true,B.true,changeT);
%         [Pmat,~] = hmmCondDecode(numStates,seq,pStates,pStates2);
%         jointm_temp = squeeze(sum(Pmat(:,:,:,2:end),4)).*permute(repmat(pStates,[1 1 l]),[2 3 1]);
%         jointm = squeeze(sum(jointm_temp(:,:,2:end),3));
%         tempvar(i) = sqrt(sum(sum(jointm)) - tempmu(i)^2); 
        % condition on z3=1   z8=1
        % Paths = Paths(Paths(:,20)==5,:);
        % allPath = allPath(allPath(:,8)==1,:); 
        %first 100 
        Probs = zeros(1,size(Paths,1));
        for j=1:size(Paths,1)
            Probs(j) = hmmpathprob_changeT_ds(Paths(j,:)+1,seqs(i,:),T.true,B.true,logpseq,len);
            dist(i,sum(Paths(j,:)~=0)+1)=dist(i,sum(Paths(j,:)~=0)+1)+Probs(j);
        end  
        ProbsAll(i,:) = Probs;
    end
    tempM = 0:l;
    dist = dist./repmat(sum(dist,2),1,size(dist,2));
    realmu = dist*tempM';
    %scatter(tempmu,realmu);
    realvar = sqrt(dist*(tempM.^2)'-realmu.^2);
    %scatter(realvar,tempvar);
    
    SortedProbsAll=sort(ProbsAll,2,'descend');
    
    [currentState, logP] = hmmviterbi_topK(seqs(1,:),wrap_T(0.35,6),B.true,5);
    
    hmmpathprob_changeT_ds(currentState(:,1),seqs(1,:),T.true,B.true,logpseq,len);
